//! Background solving. Each request runs on its own thread with its own
//! cancellation flag. Results carry the job id, so results from superseded
//! jobs are simply ignored by the UI.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{Receiver, Sender, channel};
use std::time::{Duration, Instant};

use stellar_core::{Physics, SolveControl, SolveError, SolveEvent, StarModel, StarParams, solve_model};

/// Masses used for the computed ZAMS line in the HR diagram [M☉].
pub const ZAMS_MASSES: &[f64] = &[
    1.0, 0.85, 0.7, 0.6, 0.5, 0.4, 0.3, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 13.0, 16.0, 20.0, 25.0,
    30.0, 40.0,
];

/// Composition and convection parameters a ZAMS line depends on.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZamsKey {
    pub x: f64,
    pub z: f64,
    pub alpha: f64,
    pub f_cno: f64,
}

impl ZamsKey {
    pub fn of(p: &StarParams) -> Self {
        Self { x: p.x, z: p.z, alpha: p.alpha_mlt, f_cno: p.f_cno }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ZamsPoint {
    pub mass: f64,
    pub teff: f64,
    pub luminosity: f64,
}

pub enum Msg {
    Progress { job: u64, event: SolveEvent },
    Finished { job: u64, result: Box<Result<StarModel, SolveError>>, elapsed: Duration },
    Zams { key: ZamsKey, point: ZamsPoint },
    ZamsDone { key: ZamsKey },
}

struct Ctl {
    job: u64,
    cancel: Arc<AtomicBool>,
    tx: Sender<Msg>,
    ctx: egui::Context,
    last: Option<Instant>,
}

impl SolveControl for Ctl {
    fn cancelled(&self) -> bool {
        self.cancel.load(Ordering::Relaxed)
    }

    fn report(&mut self, event: SolveEvent) {
        // Stage changes always go through; iterations are throttled.
        let is_stage = matches!(event, SolveEvent::Stage(_));
        let due = self.last.is_none_or(|t| t.elapsed() > Duration::from_millis(30));
        if is_stage || due {
            self.last = Some(Instant::now());
            let _ = self.tx.send(Msg::Progress { job: self.job, event });
            self.ctx.request_repaint();
        }
    }
}

pub struct Worker {
    physics: Arc<Physics>,
    ctx: egui::Context,
    tx: Sender<Msg>,
    rx: Receiver<Msg>,
    next_job: u64,
    current: Option<(u64, Arc<AtomicBool>)>,
    zams_cancel: Option<Arc<AtomicBool>>,
}

impl Worker {
    pub fn new(ctx: egui::Context, physics: Arc<Physics>) -> Self {
        let (tx, rx) = channel();
        Self { physics, ctx, tx, rx, next_job: 1, current: None, zams_cancel: None }
    }

    pub fn physics(&self) -> &Arc<Physics> {
        &self.physics
    }

    /// Start solving `params`. This cancels any job still running.
    pub fn submit(&mut self, params: StarParams, warm: Option<Arc<StarModel>>) -> u64 {
        self.cancel();
        let job = self.next_job;
        self.next_job += 1;
        let cancel = Arc::new(AtomicBool::new(false));
        self.current = Some((job, cancel.clone()));
        let (tx, ctx, physics) = (self.tx.clone(), self.ctx.clone(), self.physics.clone());
        std::thread::spawn(move || {
            let t0 = Instant::now();
            let mut ctl = Ctl { job, cancel, tx: tx.clone(), ctx: ctx.clone(), last: None };
            let result = solve_model(&params, &physics, warm.as_deref(), &mut ctl);
            let _ = tx.send(Msg::Finished { job, result: Box::new(result), elapsed: t0.elapsed() });
            ctx.request_repaint();
        });
        job
    }

    pub fn cancel(&mut self) {
        if let Some((_, flag)) = self.current.take() {
            flag.store(true, Ordering::Relaxed);
        }
    }

    pub fn is_current(&self, job: u64) -> bool {
        self.current.as_ref().is_some_and(|(j, _)| *j == job)
    }

    pub fn finish(&mut self, job: u64) {
        if self.is_current(job) {
            self.current = None;
        }
    }

    pub fn busy(&self) -> bool {
        self.current.is_some()
    }

    /// Compute a ZAMS line for the composition in `base`, warm-starting from
    /// `seed` (a converged model with the same composition).
    pub fn submit_zams(&mut self, base: StarParams, seed: Arc<StarModel>) {
        if let Some(flag) = self.zams_cancel.take() {
            flag.store(true, Ordering::Relaxed);
        }
        let cancel = Arc::new(AtomicBool::new(false));
        self.zams_cancel = Some(cancel.clone());
        let (tx, ctx, physics) = (self.tx.clone(), self.ctx.clone(), self.physics.clone());
        let key = ZamsKey::of(&base);
        std::thread::spawn(move || {
            let mut ctl = Ctl {
                job: 0,
                cancel: cancel.clone(),
                tx: tx.clone(),
                ctx: ctx.clone(),
                last: Some(Instant::now() + Duration::from_secs(3600)),
            };
            // Walk down from 1 M☉, then restart from 1 M☉ and walk up, so
            // every step is a small continuation step.
            let mut prev: Option<StarModel> = Some((*seed).clone());
            for (i, &mass) in ZAMS_MASSES.iter().enumerate() {
                if cancel.load(Ordering::Relaxed) {
                    return;
                }
                if i > 0 && mass > ZAMS_MASSES[i - 1] && ZAMS_MASSES[i - 1] < 1.0 {
                    prev = Some((*seed).clone());
                }
                let p = StarParams { mass, n_points: 400, ..base.clone() };
                match solve_model(&p, &physics, prev.as_ref(), &mut ctl) {
                    Ok(m) => {
                        let point = ZamsPoint { mass, teff: m.globals.teff, luminosity: m.globals.luminosity };
                        let _ = tx.send(Msg::Zams { key, point });
                        ctx.request_repaint();
                        prev = Some(m);
                    }
                    Err(SolveError::Cancelled) => return,
                    Err(_) => {}
                }
            }
            let _ = tx.send(Msg::ZamsDone { key });
            ctx.request_repaint();
        });
    }

    pub fn poll(&self) -> Vec<Msg> {
        self.rx.try_iter().collect()
    }
}
