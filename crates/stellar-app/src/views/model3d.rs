//! 3D cut-away model of the star, with animated energy transport.
//!
//! * The star is a sphere with one octant removed (x, y, z > 0). The three
//!   exposed quarter-discs are coloured by the radial profile (T or ρ).
//! * The energy-generating core (or the convective core) is drawn as an uncut
//!   ball sitting in the opened corner.
//! * Radiative zones: photon packets random-walk outward across the cut faces.
//! * Convective zones: circulating cells, with rising hot gas (light) and
//!   sinking cool gas (blue). A convective envelope shows granulation on the
//!   surface, and a convective core "boils".
//!
//! Rendering is done on the CPU: triangles are rotated, shaded, depth-sorted
//! (painter's algorithm) and handed to egui as one mesh. Overlays (photons,
//! cells, zone boundaries) use an exact visibility test against the cut-away
//! geometry. The projection is orthographic, so zooming in on the thin outer
//! layers keeps them undistorted.

use std::f32::consts::{FRAC_PI_2, TAU};
use std::sync::Arc;

use egui::{Align2, Color32, FontId, Mesh, PointerButton, Pos2, Rect, Sense, Shape, Stroke, Vec2, pos2, vec2};
use stellar_core::StarModel;

use super::cross_section::{Field, colorbar, interp};
use super::{ViewCtx, fmt_value};
use crate::evolution::RunRanges;
use crate::theme::{self, Ramp};

type V3 = [f32; 3];

fn dot(a: V3, b: V3) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn scale(a: V3, s: f32) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn normalize(a: V3) -> V3 {
    let l = dot(a, a).sqrt().max(1e-12);
    scale(a, 1.0 / l)
}

/// Rotation R = Rx(pitch) · Ry(yaw); view = R · model. +z points to the viewer.
#[derive(Clone, Copy)]
struct Rot([[f32; 3]; 3]);

impl Rot {
    fn new(yaw: f32, pitch: f32) -> Self {
        let (sy, cy) = yaw.sin_cos();
        let (sp, cp) = pitch.sin_cos();
        Self([[cy, 0.0, sy], [sp * sy, cp, -sp * cy], [-cp * sy, sp, cp * cy]])
    }
    fn apply(&self, p: V3) -> V3 {
        [dot(self.0[0], p), dot(self.0[1], p), dot(self.0[2], p)]
    }
    /// Direction towards the viewer, in model coordinates (Rᵀ ẑ).
    fn view_dir(&self) -> V3 {
        self.0[2]
    }
}

/// The removed octant's three bounding faces: plane normal, and the two axis
/// directions spanning the quarter disc (rim from `a` to `b`).
const FACES: [(V3, V3, V3); 3] = [
    ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]), // plane z = 0
    ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]), // plane x = 0
    ([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]), // plane y = 0
];
/// Angular subdivisions per octant edge (shared by sphere patches and face rims).
const N: usize = 24;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FaceQuantity {
    Temperature,
    Density,
    /// Hydrogen mass fraction: shows the helium core growing during evolution.
    Hydrogen,
}

/// Seconds of animation time per pulsation cycle (the real period is shown in the legend).
const PULSE_VIS_PERIOD: f32 = 1.6;

/// Radial pulsation of a model inside the instability strip.
struct ScenePulse {
    /// Eigenfunction ξ = δr/r (1 at the surface) sampled on normalised radius.
    rn: Vec<f32>,
    xi: Vec<f32>,
    period_days: f64,
    kind: stellar_core::pulsation::PulsatorType,
}

impl ScenePulse {
    fn xi_at(&self, r: f32) -> f32 {
        if r <= self.rn[0] {
            return self.xi[0];
        }
        let i = self.rn.partition_point(|&x| x < r).clamp(1, self.rn.len() - 1);
        let w = ((r - self.rn[i - 1]) / (self.rn[i] - self.rn[i - 1]).max(1e-9)).clamp(0.0, 1.0);
        self.xi[i - 1] + w * (self.xi[i] - self.xi[i - 1])
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Surface,
    Face,
    Core,
}

#[derive(Clone, Copy)]
struct Vert {
    p: V3,
    n: V3,
    c: [f32; 3],
    r: f32,
    kind: Kind,
}

struct Scene {
    /// (model, colouring, core as sphere, colour range bits): rebuilt when any changes.
    key: SceneKey,
    verts: Vec<Vert>,
    tris: Vec<[u32; 3]>,
    field: Field,
    rn: Vec<f64>,
    r_core: f32,
    core_label: &'static str,
    core_convective: bool,
    envelope_convective: bool,
    /// Convective intervals [a, b] in r/R outside the core ball.
    convective: Vec<(f32, f32)>,
    /// Radiative intervals [a, b] in r/R outside the core ball.
    radiative: Vec<(f32, f32)>,
    /// Zone boundaries (r/R) to outline on the faces.
    boundaries: Vec<f32>,
    /// Fundamental radial mode, if the model lies in the instability strip.
    pulsation: Option<ScenePulse>,
}

pub struct Model3dState {
    pub yaw: f32,
    pub pitch: f32,
    pub zoom: f32,
    pub pan: Vec2,
    pub animate: bool,
    pub speed: f32,
    pub auto_rotate: bool,
    pub quantity: FaceQuantity,
    pub core_sphere: bool,
    /// Animate radial pulsation when the star is in the instability strip.
    pub pulsate: bool,
    /// Displayed relative radius amplitude of the pulsation (exaggerated).
    pub amplitude: f32,
    /// During an evolution run: draw the star to scale, relative to its largest radius in the run.
    pub true_size: bool,
    sim_time: f64,
    /// Zoom onto the envelope on the next frame (needs the canvas size).
    pending_surface_zoom: bool,
    scene: Option<Scene>,
}

const DEFAULT_YAW: f32 = -0.44;
const DEFAULT_PITCH: f32 = 0.35;

impl Default for Model3dState {
    fn default() -> Self {
        Self {
            yaw: DEFAULT_YAW,
            pitch: DEFAULT_PITCH,
            zoom: 1.0,
            pan: Vec2::ZERO,
            animate: true,
            speed: 1.0,
            auto_rotate: false,
            quantity: FaceQuantity::Temperature,
            core_sphere: true,
            pulsate: true,
            amplitude: 0.06,
            true_size: false,
            sim_time: 0.0,
            pending_surface_zoom: false,
            scene: None,
        }
    }
}

impl Model3dState {
    /// Stop the animation at a given animation time (used by the capture mode).
    pub fn freeze_at(&mut self, t: f64) {
        self.sim_time = t;
        self.animate = false;
    }

    /// Request a zoom onto the outer envelope at the cut (resolved on the next
    /// frame, once the canvas size is known).
    pub fn zoom_to_surface(&mut self) {
        self.yaw = DEFAULT_YAW;
        self.pitch = DEFAULT_PITCH;
        self.pan = Vec2::ZERO;
        self.pending_surface_zoom = true;
    }

    fn reset_view(&mut self) {
        self.yaw = DEFAULT_YAW;
        self.pitch = DEFAULT_PITCH;
        self.zoom = 1.0;
        self.pan = Vec2::ZERO;
    }
}

// ---------------------------------------------------------------------------
// Colours

/// Approximate sRGB colour of a black body at temperature `t` [K].
fn blackbody(t: f64) -> [f32; 3] {
    const TABLE: &[(f64, [u8; 3])] = &[
        (1000.0, [255, 56, 0]),
        (1500.0, [255, 109, 0]),
        (2000.0, [255, 137, 18]),
        (2500.0, [255, 161, 72]),
        (3000.0, [255, 180, 107]),
        (3500.0, [255, 196, 137]),
        (4000.0, [255, 209, 163]),
        (4500.0, [255, 219, 186]),
        (5000.0, [255, 228, 206]),
        (5500.0, [255, 236, 224]),
        (6000.0, [255, 243, 239]),
        (6500.0, [255, 249, 253]),
        (7000.0, [245, 243, 255]),
        (8000.0, [227, 233, 255]),
        (10000.0, [204, 219, 255]),
        (15000.0, [179, 204, 255]),
        (20000.0, [168, 197, 255]),
        (30000.0, [159, 191, 255]),
        (40000.0, [155, 188, 255]),
    ];
    let t = t.clamp(TABLE[0].0, TABLE[TABLE.len() - 1].0);
    let i = TABLE.partition_point(|(tt, _)| *tt < t).clamp(1, TABLE.len() - 1);
    let (t0, c0) = TABLE[i - 1];
    let (t1, c1) = TABLE[i];
    let w = ((t - t0) / (t1 - t0)) as f32;
    std::array::from_fn(|k| (c0[k] as f32 + w * (c1[k] as f32 - c0[k] as f32)) / 255.0)
}

fn to_f3(c: Color32) -> [f32; 3] {
    [c.r() as f32 / 255.0, c.g() as f32 / 255.0, c.b() as f32 / 255.0]
}

fn to_color(c: [f32; 3]) -> Color32 {
    let q = |v: f32| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
    Color32::from_rgb(q(c[0]), q(c[1]), q(c[2]))
}

/// Light from the upper left, in view space.
const LIGHT: V3 = [-0.42, 0.55, 0.72];

/// Smooth, evolving cellular pattern used for granulation and boiling.
fn granulation(p: V3, t: f32, freq: f32) -> f32 {
    let a = (freq * p[0] + 0.9 * t).sin() * (freq * p[1] - 0.7 * t + 1.3).sin() * (freq * p[2] + 0.5 * t + 2.1).sin();
    let b = (1.7 * freq * (p[0] + p[1]) - 1.1 * t).sin() * (1.7 * freq * (p[1] - p[2]) + 0.8 * t + 0.4).sin();
    a + 0.5 * b
}

// ---------------------------------------------------------------------------
// Scene construction (model space, cached per model)

type SceneKey = (usize, FaceQuantity, bool, (u64, u64));

/// Colour range for the cut faces with a fixed (run-wide) scale, in scale units.
fn fixed_face_range(quantity: FaceQuantity, fixed: Option<&RunRanges>) -> Option<(f64, f64)> {
    let r = fixed?;
    let range = match quantity {
        FaceQuantity::Temperature => r.log_t,
        FaceQuantity::Density => r.log_rho,
        FaceQuantity::Hydrogen => crate::evolution::Range { lo: 0.0, hi: r.x_max },
    };
    range.is_valid().then_some((range.lo, range.hi))
}

fn scene_key(model: &Arc<StarModel>, quantity: FaceQuantity, core_sphere: bool, fixed: Option<&RunRanges>) -> SceneKey {
    let bits = fixed_face_range(quantity, fixed).map(|(lo, hi)| (lo.to_bits(), hi.to_bits())).unwrap_or((0, 0));
    (Arc::as_ptr(model) as usize, quantity, core_sphere, bits)
}

fn build_scene(model: &Arc<StarModel>, quantity: FaceQuantity, core_sphere: bool, fixed: Option<&RunRanges>) -> Scene {
    let p = &model.profile;
    let r_surf = *p.r.last().unwrap();
    let rn: Vec<f64> = p.r.iter().map(|r| r / r_surf).collect();
    // Physical metaphor: hotter / denser = brighter, in both themes.
    let field = match quantity {
        FaceQuantity::Temperature => Field::new(&p.t, Ramp::orange(true)),
        FaceQuantity::Density => Field::new(&p.rho, Ramp::blue(true)),
        FaceQuantity::Hydrogen => {
            let x_max = p.x.iter().cloned().fold(0.0, f64::max);
            // Hydrogen-rich gas dark, the burnt (helium) core light: the core
            // stands out and the overlays keep their contrast.
            Field::linear(&p.x, Ramp::aqua(false), x_max)
        }
    };
    let field = match fixed_face_range(quantity, fixed) {
        Some((lo, hi)) => field.with_range(lo, hi),
        None => field,
    };
    let g = &model.globals;

    // The ball: the convective core if there is one (and it is not the whole
    // star), otherwise the region producing 90% of the luminosity.
    let l_surf = p.l[p.len() - 1];
    let l_norm: Vec<f64> = p.l.iter().map(|l| l / l_surf).collect();
    let r90 = {
        let i = l_norm.partition_point(|&l| l < 0.9).clamp(1, p.len() - 1);
        let w = (0.9 - l_norm[i - 1]) / (l_norm[i] - l_norm[i - 1]).max(1e-30);
        rn[i - 1] + w * (rn[i] - rn[i - 1])
    };
    let (r_core, core_label) = if !core_sphere {
        (0.0, "")
    } else if g.conv_core_q > 0.0 && g.conv_core_r < 0.95 {
        (g.conv_core_r, "convective core")
    } else {
        (r90, "core producing 90% of L")
    };
    let r_core = r_core as f32;
    let core_convective = r_core > 0.0 && {
        let k = rn.partition_point(|&r| r < 0.5 * r_core as f64).min(p.len() - 1);
        p.convective[k]
    };

    // Zones outside the ball.
    let mut convective = Vec::new();
    for z in &model.zones {
        let (a, b) = ((z.r_inner as f32).max(r_core), z.r_outer as f32);
        if b > a + 1e-6 {
            convective.push((a, b));
        }
    }
    let mut radiative = Vec::new();
    let mut cursor = r_core;
    for &(a, b) in &convective {
        if a > cursor + 1e-6 {
            radiative.push((cursor, a));
        }
        cursor = cursor.max(b);
    }
    if cursor < 1.0 - 1e-6 {
        radiative.push((cursor, 1.0));
    }
    let envelope_convective = convective.last().is_some_and(|&(_, b)| b > 0.999);
    let mut boundaries: Vec<f32> = model
        .zones
        .iter()
        .flat_map(|z| [z.r_inner as f32, z.r_outer as f32])
        .filter(|&r| r > r_core + 1e-4 && r < 0.9995)
        .collect();
    boundaries.dedup();

    let mut verts = Vec::new();
    let mut tris = Vec::new();

    // Outer surface: the seven octants that are not cut away.
    let surface_color = blackbody(g.teff);
    for sx in [-1.0_f32, 1.0] {
        for sy in [-1.0_f32, 1.0] {
            for sz in [-1.0_f32, 1.0] {
                if sx > 0.0 && sy > 0.0 && sz > 0.0 {
                    continue;
                }
                octant_patch(&mut verts, &mut tris, [sx, sy, sz], 1.0, surface_color, Kind::Surface);
            }
        }
    }

    // Cut faces: quarter discs from the ball (or centre) to the surface.
    let radii = ring_radii(&rn, r_core as f64, &boundaries);
    for (normal, a, b) in FACES {
        let base = verts.len() as u32;
        for &r in &radii {
            let c = to_f3(field.color_at(&rn, r));
            for j in 0..=N {
                let s = j as f32 / N as f32;
                let d = normalize(add(scale(a, 1.0 - s), scale(b, s)));
                verts.push(Vert { p: scale(d, r as f32), n: normal, c, r: r as f32, kind: Kind::Face });
            }
        }
        let row = (N + 1) as u32;
        for i in 0..(radii.len() - 1) as u32 {
            for j in 0..N as u32 {
                let v00 = base + i * row + j;
                tris.push([v00, v00 + row, v00 + row + 1]);
                tris.push([v00, v00 + row + 1, v00 + 1]);
            }
        }
    }

    // The core ball: only its octant inside the opened corner is ever visible.
    if r_core > 0.0 {
        let c = to_f3(field.color_at(&rn, 0.5 * r_core as f64));
        octant_patch(&mut verts, &mut tris, [1.0, 1.0, 1.0], r_core, c, Kind::Core);
    }

    let pulsation = model.pulsation.as_ref().and_then(|mode| {
        let kind = stellar_core::pulsation::classify(g.teff, g.luminosity)?;
        let xi_s = *mode.xi.last()?;
        (mode.xi.len() == rn.len()).then(|| ScenePulse {
            rn: rn.iter().map(|&r| r as f32).collect(),
            xi: mode.xi.iter().map(|&x| (x / xi_s) as f32).collect(),
            period_days: mode.period_days(),
            kind,
        })
    });

    Scene {
        key: scene_key(model, quantity, core_sphere, fixed),
        verts,
        tris,
        field,
        rn,
        r_core,
        core_label,
        core_convective,
        envelope_convective,
        convective,
        radiative,
        boundaries,
        pulsation,
    }
}

/// One octant of a sphere of radius `radius`, from a subdivided octahedron face.
fn octant_patch(verts: &mut Vec<Vert>, tris: &mut Vec<[u32; 3]>, s: V3, radius: f32, c: [f32; 3], kind: Kind) {
    let a = [s[0], 0.0, 0.0];
    let b = [0.0, s[1], 0.0];
    let cc = [0.0, 0.0, s[2]];
    let base = verts.len() as u32;
    // Row i holds N − i + 1 vertices, so row i starts after Σ_{k<i} (N − k + 1).
    let index = |i: usize, j: usize| -> u32 { base + (i * (N + 1) - i * (i.saturating_sub(1)) / 2 + j) as u32 };
    for i in 0..=N {
        for j in 0..=(N - i) {
            let k = N - i - j;
            let p = add(add(scale(a, i as f32), scale(b, j as f32)), scale(cc, k as f32));
            let n = normalize(p);
            debug_assert_eq!(index(i, j), verts.len() as u32);
            verts.push(Vert { p: scale(n, radius), n, c, r: radius, kind });
        }
    }
    for i in 0..N {
        for j in 0..(N - i) {
            tris.push([index(i, j), index(i + 1, j), index(i, j + 1)]);
            if i + j + 1 < N {
                tris.push([index(i + 1, j), index(i + 1, j + 1), index(i, j + 1)]);
            }
        }
    }
}

/// Ring radii for the cut faces: uniform steps, refined towards the surface
/// (using the solver mesh, which is dense there) and at zone boundaries.
fn ring_radii(rn: &[f64], r0: f64, boundaries: &[f32]) -> Vec<f64> {
    const UNIFORM: usize = 120;
    let mut v: Vec<f64> = (0..=UNIFORM).map(|i| r0 + (1.0 - r0) * i as f64 / UNIFORM as f64).collect();
    v.extend(rn.iter().copied().filter(|&r| r > r0 && r > 0.85));
    v.extend(boundaries.iter().map(|&b| b as f64));
    v.sort_by(f64::total_cmp);
    let mut out: Vec<f64> = Vec::with_capacity(v.len());
    for r in v {
        // Keep spacing relative to the remaining depth, so the outermost
        // layers stay resolved when zooming in.
        let min_gap = (2e-3 * (1.0 - r)).max(2e-5);
        if out.last().is_none_or(|&last| r - last > min_gap) {
            out.push(r);
        }
    }
    if *out.last().unwrap() < 1.0 {
        out.push(1.0);
    } else {
        *out.last_mut().unwrap() = 1.0;
    }
    out
}

// ---------------------------------------------------------------------------
// Visibility of a point on a cut face (orthographic view along `d`).

fn visible(p: V3, d: V3, r_core: f32) -> bool {
    // Leaving the open octant through one of its planes before the sphere
    // surface means the ray enters solid material.
    let pd = dot(p, d);
    let t_exit = -pd + (pd * pd - (dot(p, p) - 1.0)).max(0.0).sqrt();
    for i in 0..3 {
        if d[i] < 0.0 {
            if p[i] <= 1e-6 {
                return false;
            }
            if -p[i] / d[i] < t_exit {
                return false;
            }
        }
    }
    if r_core > 0.0 {
        let c = dot(p, p) - r_core * r_core;
        let disc = pd * pd - c;
        if disc > 0.0 {
            let t1 = -pd - disc.sqrt();
            if t1 > 1e-4 && t1 < t_exit {
                return false;
            }
        }
    }
    true
}

/// Point on face `f` at radius `r` and angle `phi` ∈ [0, π/2].
fn face_point(f: usize, r: f32, phi: f32) -> V3 {
    let (_, a, b) = FACES[f];
    scale(add(scale(a, phi.cos()), scale(b, phi.sin())), r)
}

/// Deterministic pseudo-random number in [0, 1) from an integer seed.
fn hash(i: u32, salt: u32) -> f32 {
    let mut x = i.wrapping_mul(0x9E37_79B9) ^ salt.wrapping_mul(0x85EB_CA6B);
    x ^= x >> 15;
    x = x.wrapping_mul(0x2C1B_3C6D);
    x ^= x >> 12;
    x = x.wrapping_mul(0x297A_2D39);
    x ^= x >> 15;
    (x >> 8) as f32 / (1u32 << 24) as f32
}

fn smoothstep(e0: f32, e1: f32, x: f32) -> f32 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Photons: small white sparks with a soft glow and a trail.
const PHOTON: Color32 = Color32::from_rgb(0xff, 0xff, 0xf4);
/// Convective parcels: outlined discs, warm amber rising and blue sinking.
const RISING: Color32 = Color32::from_rgb(0xff, 0xc8, 0x4a);
const SINKING: Color32 = Color32::from_rgb(0x3d, 0x8b, 0xff);

/// A photon spark: bright core with a faint halo.
fn photon(painter: &egui::Painter, p: Pos2, r: f32, alpha: f32) {
    painter.circle_filled(p, r * 1.9, PHOTON.gamma_multiply(0.22 * alpha));
    painter.circle_filled(p, r * 0.75, PHOTON.gamma_multiply(alpha));
}

// ---------------------------------------------------------------------------

pub fn show(ui: &mut egui::Ui, v: &ViewCtx<'_>, state: &mut Model3dState) {
    let dark = v.dark;

    let key = scene_key(v.model, state.quantity, state.core_sphere, v.fixed);
    if state.scene.as_ref().is_none_or(|s| s.key != key) {
        state.scene = Some(build_scene(v.model, state.quantity, state.core_sphere, v.fixed));
    }
    let pulsating = state.scene.as_ref().is_some_and(|s| s.pulsation.is_some());

    // Controls, in short rows so they fit narrow windows.
    ui.horizontal(|ui| {
        ui.label("Cut faces:");
        ui.radio_value(&mut state.quantity, FaceQuantity::Temperature, "temperature");
        ui.radio_value(&mut state.quantity, FaceQuantity::Density, "density");
        ui.radio_value(&mut state.quantity, FaceQuantity::Hydrogen, "hydrogen X")
            .on_hover_text("hydrogen mass fraction: the helium core grows as the star evolves");
        ui.separator();
        ui.checkbox(&mut state.core_sphere, "core as sphere");
        ui.separator();
        ui.checkbox(&mut state.animate, "animate");
        ui.add(egui::Slider::new(&mut state.speed, 0.1..=5.0).logarithmic(true).text("speed").max_decimals(1));
        ui.checkbox(&mut state.auto_rotate, "auto-rotate");
        ui.separator();
        ui.checkbox(&mut state.true_size, "true size").on_hover_text(
            "during an evolution run: draw the star to scale, relative to its largest radius in the run \
             (zoom in to see small stages)",
        );
        ui.add_enabled(pulsating, egui::Checkbox::new(&mut state.pulsate, "pulsation"))
            .on_hover_text("radial pulsation, shown when the star lies in the instability strip")
            .on_disabled_hover_text("this model is outside the instability strip: no classical pulsation");
    });
    ui.horizontal(|ui| {
        ui.add(egui::Slider::new(&mut state.zoom, 0.3..=400.0).logarithmic(true).text("zoom").max_decimals(1));
        if ui.button("Zoom to surface layers").on_hover_text("magnify the outer envelope at the cut").clicked() {
            state.zoom_to_surface();
        }
        if ui.button("Reset view").clicked() {
            state.reset_view();
        }
        ui.label(
            egui::RichText::new("drag: rotate · scroll/pinch: zoom · right- or shift-drag: pan · double-click: reset")
                .small()
                .color(theme::muted(dark)),
        );
    });

    let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
    let rect = resp.rect;
    let painter = painter.with_clip_rect(rect);
    let base_px = 0.40 * rect.width().min(rect.height());
    let center = rect.center();

    // Interaction.
    let (dt, shift) = ui.input(|i| (i.stable_dt.min(0.1), i.modifiers.shift));
    let d = resp.drag_delta();
    if resp.dragged_by(PointerButton::Primary) && !shift {
        state.yaw += d.x * 0.008;
        state.pitch = (state.pitch + d.y * 0.008).clamp(-1.55, 1.55);
    } else if resp.dragged_by(PointerButton::Secondary)
        || resp.dragged_by(PointerButton::Middle)
        || (resp.dragged_by(PointerButton::Primary) && shift)
    {
        state.pan += d;
    }
    if resp.hovered() {
        let (scroll, pinch) = ui.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
        let factor = pinch * (scroll * 0.003).exp();
        if (factor - 1.0).abs() > 1e-4 {
            let new_zoom = (state.zoom * factor).clamp(0.3, 400.0);
            let f = new_zoom / state.zoom;
            if let Some(m) = resp.hover_pos() {
                state.pan = (m - center) - f * ((m - center) - state.pan);
            }
            state.zoom = new_zoom;
        }
    }
    if resp.double_clicked() {
        state.reset_view();
    }
    if state.auto_rotate {
        state.yaw += dt * 0.25;
    }
    if state.animate {
        state.sim_time += (dt * state.speed) as f64;
    }

    let scene = state.scene.as_ref().unwrap();
    // True size: the star is drawn to scale relative to its largest radius in
    // the run (with a faint outline of that size).
    let (size, r_max) = match v.run.filter(|_| state.true_size) {
        Some(rr) if rr.radius.is_valid() && rr.radius.hi > 0.0 => {
            ((v.model.globals.radius / rr.radius.hi).clamp(1e-3, 1.0) as f32, rr.radius.hi)
        }
        _ => (1.0, v.model.globals.radius),
    };
    let rot = Rot::new(state.yaw, state.pitch);

    // "Zoom to surface layers": centre on the middle of the envelope at the cut.
    if std::mem::take(&mut state.pending_surface_zoom) {
        let base = scene.convective.last().filter(|z| z.1 > 0.999).map(|z| z.0).unwrap_or(0.95);
        let depth = (1.0 - base).max(0.004);
        state.zoom = (0.45 / depth).clamp(1.0, 400.0);
        let target = face_point(0, 0.5 * (base + 1.0), std::f32::consts::FRAC_PI_4);
        let q = rot.apply(target);
        state.pan = -vec2(q[0], -q[1]) * base_px * state.zoom * size;
    }

    let scale_px = base_px * state.zoom * size;
    let origin = center + state.pan;
    let project = |p: V3| -> (Pos2, f32) {
        let q = rot.apply(p);
        (pos2(origin.x + q[0] * scale_px, origin.y - q[1] * scale_px), q[2])
    };
    let t = state.sim_time as f32;
    // Radial pulsation: each layer moves by ξ(r) (1 at the surface) times the
    // displayed amplitude; layers are shown at r · (1 + A ξ(r) sin φ).
    let pulse = scene.pulsation.as_ref().filter(|_| state.pulsate);
    let pulse_s = (TAU * t / PULSE_VIS_PERIOD).sin();
    let amp = state.amplitude;
    let disp = |r: f32| -> f32 { pulse.map_or(1.0, |p| 1.0 + amp * p.xi_at(r) * pulse_s) };
    // Project a point of the static geometry at its pulsating position.
    let proj = |p: V3| project(scale(p, disp(dot(p, p).sqrt())));

    if size < 0.995 {
        let full = base_px * state.zoom;
        painter.circle_stroke(origin, full, Stroke::new(1.0, theme::muted(dark).gamma_multiply(0.6)));
    }

    // Soft ground shadow (light theme only; invisible on dark backgrounds).
    if !dark && state.pan.length() < 1.0 && state.zoom < 1.5 {
        let c = origin + vec2(0.0, 1.12 * scale_px);
        let (rx, ry) = (0.85 * scale_px, 0.10 * scale_px);
        let mut mesh = Mesh::default();
        mesh.colored_vertex(c, Color32::from_black_alpha(55));
        for k in 0..=48 {
            let a = TAU * k as f32 / 48.0;
            mesh.colored_vertex(c + vec2(a.cos() * rx, a.sin() * ry), Color32::TRANSPARENT);
        }
        for k in 1..=48u32 {
            mesh.add_triangle(0, k, k + 1);
        }
        painter.add(Shape::mesh(mesh));
    }

    // Transform, shade and depth-sort the triangles.
    let light = normalize(LIGHT);
    let half = normalize(add(light, [0.0, 0.0, 1.0]));
    let mut screen = Vec::with_capacity(scene.verts.len());
    let mut depth = Vec::with_capacity(scene.verts.len());
    let mut nz = Vec::with_capacity(scene.verts.len());
    let mut colors = Vec::with_capacity(scene.verts.len());
    // The surface is hottest and brightest near maximum compression.
    let pulse_light = if pulse.is_some() { 1.0 - 2.5 * amp * pulse_s } else { 1.0 };
    for vtx in &scene.verts {
        let (s, z) = project(scale(vtx.p, disp(vtx.r)));
        let n = rot.apply(vtx.n);
        screen.push(s);
        depth.push(z);
        nz.push(n[2]);
        let diffuse = dot(n, light).max(0.0);
        let spec = dot(n, half).max(0.0).powi(40);
        let c = match vtx.kind {
            Kind::Surface => {
                let gran = if scene.envelope_convective { 1.0 + 0.07 * granulation(vtx.p, t * 0.6, 18.0) } else { 1.0 };
                let k = (0.38 + 0.62 * diffuse) * gran * pulse_light;
                vtx.c.map(|c| c * k + 0.30 * spec)
            }
            Kind::Face => {
                let k = 0.80 + 0.20 * diffuse;
                vtx.c.map(|c| c * k)
            }
            Kind::Core => {
                let boil = if scene.core_convective { 1.0 + 0.12 * granulation(vtx.p, t * 1.2, 30.0) } else { 1.0 };
                let k = (0.45 + 0.55 * diffuse) * boil;
                vtx.c.map(|c| c * k + 0.35 * spec)
            }
        };
        colors.push(to_color(c));
    }
    let mut order: Vec<(f32, u32)> = scene
        .tris
        .iter()
        .enumerate()
        .filter(|(_, tri)| tri.iter().map(|&i| nz[i as usize]).sum::<f32>() > -0.05)
        .map(|(k, tri)| (tri.iter().map(|&i| depth[i as usize]).sum::<f32>(), k as u32))
        .collect();
    order.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let mut mesh = Mesh::default();
    mesh.reserve_vertices(order.len() * 3);
    mesh.reserve_triangles(order.len());
    for (n, &(_, k)) in order.iter().enumerate() {
        for &i in &scene.tris[k as usize] {
            mesh.colored_vertex(screen[i as usize], colors[i as usize]);
        }
        let b = 3 * n as u32;
        mesh.add_triangle(b, b + 1, b + 2);
    }
    painter.add(Shape::mesh(mesh));

    // Overlays on the cut faces.
    let view_dir = rot.view_dir();
    let r_core = scene.r_core;
    let vis = |p: V3| visible(p, view_dir, r_core);
    let dot_r = (2.2 * state.zoom.sqrt()).clamp(2.0, 6.0);

    // Zone boundaries.
    let ink = Color32::from_white_alpha(150);
    for &rb in &scene.boundaries {
        for f in 0..3 {
            let pts: Vec<(Pos2, bool)> = (0..=48)
                .map(|k| {
                    let p = face_point(f, rb, FRAC_PI_2 * k as f32 / 48.0);
                    (proj(p).0, vis(p))
                })
                .collect();
            for w in pts.windows(2) {
                if w[0].1 && w[1].1 {
                    painter.line_segment([w[0].0, w[1].0], Stroke::new(1.2, ink));
                }
            }
        }
    }

    // Photons random-walking outward through radiative zones.
    let total_area: f32 = scene.radiative.iter().map(|&(a, b)| b * b - a * a).sum();
    let n_photons = 180usize;
    let mut idx = 0u32;
    for &(a, b) in &scene.radiative {
        let count = ((b * b - a * a) / total_area.max(1e-6) * n_photons as f32).round() as u32;
        for _ in 0..count.max(4) {
            idx += 1;
            let f = (idx % 3) as usize;
            let speed = 0.03 + 0.04 * hash(idx, 1);
            let phase = hash(idx, 2);
            let phi0 = 0.08 + (FRAC_PI_2 - 0.16) * hash(idx, 3);
            let walk = |tt: f32| -> V3 {
                let frac = (phase + tt * speed).fract();
                let r = a + (b - a) * frac + 0.006 * (9.3 * tt + 6.0 * hash(idx, 4)).sin();
                let jitter =
                    0.016 * ((5.1 * tt + 6.0 * hash(idx, 5)).sin() + 0.5 * (11.7 * tt + 6.0 * hash(idx, 6)).sin());
                let phi = (phi0 + jitter / r.max(0.05)).clamp(0.02, FRAC_PI_2 - 0.02);
                face_point(f, r.clamp(a, b), phi)
            };
            let frac = (phase + t * speed).fract();
            let alpha = smoothstep(0.0, 0.06, frac) * smoothstep(1.0, 0.94, frac);
            let p_now = walk(t);
            if alpha <= 0.0 || !vis(p_now) {
                continue;
            }
            let (s_now, _) = proj(p_now);
            // Short trail showing the zig-zag path.
            let mut prev = s_now;
            for k in 1..=5 {
                let pk = walk(t - 0.12 * k as f32);
                if !vis(pk) || (phase + (t - 0.12 * k as f32) * speed).fract() > frac {
                    break;
                }
                let sk = proj(pk).0;
                let a_k = alpha * (1.0 - k as f32 / 6.0) * 0.6;
                painter.line_segment([prev, sk], Stroke::new(1.0, PHOTON.gamma_multiply(a_k)));
                prev = sk;
            }
            photon(&painter, s_now, dot_r, alpha);
        }
    }

    // Convection cells.
    let omega = 1.6;
    let mut cell_id = 0u32;
    for &(a, b) in &scene.convective {
        let h = b - a;
        let rows = ((h / 0.16).ceil() as usize).max(1);
        let hr = h / rows as f32;
        for row in 0..rows {
            let rm = a + hr * (row as f32 + 0.5);
            let count = (((rm * FRAC_PI_2) / (hr * 1.15)).round() as usize).clamp(1, 40);
            let dphi = FRAC_PI_2 / count as f32;
            for f in 0..3 {
                for col in 0..count {
                    cell_id += 1;
                    let phic = dphi * (col as f32 + 0.5);
                    let dir = if (row + col + f) % 2 == 0 { 1.0 } else { -1.0 };
                    let ar = 0.40 * hr;
                    let at = 0.40 * rm * dphi;
                    let at_angle = |th: f32| -> V3 {
                        let r = rm + ar * th.cos();
                        let phi = phic + at * th.sin() / rm.max(1e-4);
                        face_point(f, r, phi)
                    };
                    // Loop outline.
                    let pts: Vec<(Pos2, bool)> = (0..=24)
                        .map(|k| {
                            let p = at_angle(TAU * k as f32 / 24.0);
                            (proj(p).0, vis(p))
                        })
                        .collect();
                    if pts.iter().all(|p| !p.1) {
                        continue;
                    }
                    for w in pts.windows(2) {
                        if w[0].1 && w[1].1 {
                            painter.line_segment([w[0].0, w[1].0], Stroke::new(1.0, Color32::from_white_alpha(70)));
                        }
                    }
                    // Parcels circulating around the loop.
                    let offset = TAU * hash(cell_id, 7);
                    for k in 0..3 {
                        let th = offset + dir * omega * t + TAU * k as f32 / 3.0;
                        let p = at_angle(th);
                        if !vis(p) {
                            continue;
                        }
                        // dr/dt ∝ −sin θ · dir: positive means the parcel rises.
                        let rising = -th.sin() * dir > 0.0;
                        let col_ = if rising { RISING } else { SINKING };
                        let sp = proj(p).0;
                        painter.circle(sp, dot_r, col_, Stroke::new(0.8, Color32::from_black_alpha(120)));
                    }
                }
            }
        }
    }

    // Colour bar and legend.
    let bar = Rect::from_min_size(
        pos2(rect.right() - 70.0, rect.top() + 40.0),
        vec2(14.0, (rect.height() * 0.45).max(120.0)),
    );
    let (title, unit) = match state.quantity {
        FaceQuantity::Temperature => ("T", "K"),
        FaceQuantity::Density => ("ρ", "g cm⁻³"),
        FaceQuantity::Hydrogen => ("X", ""),
    };
    painter.rect_filled(
        Rect::from_min_max(pos2(bar.left() - 52.0, bar.top() - 30.0), pos2(bar.right() + 10.0, bar.bottom() + 10.0)),
        6.0,
        panel_backdrop(dark),
    );
    colorbar(&painter, bar, &scene.field, title, unit, true, dark);
    let size_note = (size < 0.995).then(|| {
        format!(
            "true size: R = {} R☉; faint outline = largest radius in run ({} R☉)",
            fmt_value(v.model.globals.radius),
            fmt_value(r_max)
        )
    });
    draw_legend(&painter, rect, scene, v.model, dark, dot_r, LegendExtras { pulsing: pulse.is_some(), size_note });

    // Hover read-out: front-most triangle under the pointer.
    if let Some(m) = resp.hover_pos()
        && !resp.dragged()
    {
        let hit = order.iter().rev().find_map(|&(_, k)| {
            let [i0, i1, i2] = scene.tris[k as usize];
            let (a, b, c) = (screen[i0 as usize], screen[i1 as usize], screen[i2 as usize]);
            let den = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
            if den.abs() < 1e-9 {
                return None;
            }
            let w0 = ((b.y - c.y) * (m.x - c.x) + (c.x - b.x) * (m.y - c.y)) / den;
            let w1 = ((c.y - a.y) * (m.x - c.x) + (a.x - c.x) * (m.y - c.y)) / den;
            let w2 = 1.0 - w0 - w1;
            (w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0).then(|| {
                let vs = &scene.verts;
                let r = w0 * vs[i0 as usize].r + w1 * vs[i1 as usize].r + w2 * vs[i2 as usize].r;
                (vs[i0 as usize].kind, r)
            })
        });
        if let Some((kind, r)) = hit {
            let pr = &v.model.profile;
            let rn = &scene.rn;
            let s = r as f64;
            let log_at =
                |ys: &[f64]| 10f64.powf(interp(rn, &ys.iter().map(|y| y.max(1e-300).log10()).collect::<Vec<_>>(), s));
            let text = match kind {
                Kind::Surface => format!(
                    "Photosphere (r/R = 1)\nT_eff = {:.0} K\nsurface colour: black body at T_eff{}",
                    v.model.globals.teff,
                    if scene.envelope_convective { "\nmottling: granulation (convective envelope)" } else { "" }
                ),
                _ => {
                    let k = rn.partition_point(|&x| x < s).min(rn.len() - 1);
                    let zone = if pr.convective[k] { "convective" } else { "radiative" };
                    let head =
                        if kind == Kind::Core { format!("{} boundary\n", scene.core_label) } else { String::new() };
                    format!(
                        "{head}r/R = {:.4}\nm/M = {:.4}\nT = {} K\nρ = {} g cm⁻³\nP = {} dyn cm⁻²\nL/L_s = {:.4}\nX = {:.4}\n{zone} transport",
                        s,
                        interp(rn, &pr.q, s),
                        fmt_value(log_at(&pr.t)),
                        fmt_value(log_at(&pr.rho)),
                        fmt_value(log_at(&pr.p)),
                        interp(rn, &pr.l, s) / pr.l[pr.len() - 1],
                        interp(rn, &pr.x, s),
                    )
                }
            };
            resp.on_hover_text_at_pointer(text);
        }
    }

    if state.animate || state.auto_rotate {
        ui.ctx().request_repaint();
    }
}

/// Translucent panel colour behind overlays drawn on top of the star.
fn panel_backdrop(dark: bool) -> Color32 {
    if dark {
        Color32::from_rgba_unmultiplied(26, 26, 25, 205)
    } else {
        Color32::from_rgba_unmultiplied(252, 252, 251, 215)
    }
}

/// View-state lines for the 3D legend.
struct LegendExtras {
    /// The pulsation animation is running.
    pulsing: bool,
    /// Explanation of the true-size view, when active.
    size_note: Option<String>,
}

fn draw_legend(
    painter: &egui::Painter,
    rect: Rect,
    scene: &Scene,
    model: &StarModel,
    dark: bool,
    dot_r: f32,
    extras: LegendExtras,
) {
    let LegendExtras { pulsing, size_note } = extras;
    // Icons sit on a chip of face colour, as they appear on the star.
    let chip = scene.field.color_at(&scene.rn, 0.85);
    let ink = theme::primary(dark);
    let secondary = theme::secondary(dark);
    let font = FontId::proportional(13.0);
    let mut y = rect.bottom() - 12.0;
    let x = rect.left() + 14.0;
    let mut lines: Vec<(Option<(Color32, bool)>, String)> = Vec::new();
    if !scene.radiative.is_empty() {
        lines.push((Some((PHOTON, false)), "photon random walk: radiative diffusion outward".into()));
    }
    if !scene.convective.is_empty() {
        lines.push((Some((RISING, true)), "convection: hot gas rising".into()));
        lines.push((Some((SINKING, true)), "convection: cool gas sinking".into()));
    }
    if scene.r_core > 0.0 {
        let boil = if scene.core_convective { ", boiling = convective" } else { "" };
        lines.push((None, format!("core sphere: {} (r/R < {:.3}){boil}", scene.core_label, scene.r_core)));
    }
    let g = &model.globals;
    lines.push((
        None,
        format!(
            "surface: black-body colour at T_eff = {:.0} K{}",
            g.teff,
            if scene.envelope_convective { ", mottling = granulation" } else { "" }
        ),
    ));
    lines.push((None, "white arcs: boundaries between radiative and convective zones".into()));
    lines.push((
        None,
        "photons, convection cells and granulation are illustrative: a 1D model, not a 3D simulation".into(),
    ));
    if let Some(note) = size_note {
        lines.push((None, note));
    }
    if let Some(p) = &scene.pulsation {
        let state = if pulsing { "" } else { " (animation off)" };
        lines.push((
            None,
            format!(
                "{}: fundamental radial mode, P = {}; shown as 1 cycle per {PULSE_VIS_PERIOD} s, amplitude exaggerated{state}",
                p.kind.label(),
                crate::evolution::fmt_period(p.period_days)
            ),
        ));
    }
    let width = lines.iter().map(|(_, t)| t.chars().count()).max().unwrap_or(0) as f32 * 6.6 + 34.0;
    let height = lines.len() as f32 * 18.0 + 10.0;
    let backdrop = Rect::from_min_max(pos2(x - 6.0, y - height + 2.0), pos2(x + width, y + 6.0));
    painter.rect_filled(backdrop, 6.0, panel_backdrop(dark));
    for (icon, text) in lines.iter().rev() {
        let tx = x + 24.0;
        if let Some((c, parcel)) = icon {
            let p = pos2(x + 8.0, y - 8.0);
            painter.rect_filled(Rect::from_center_size(p - vec2(2.0, 0.0), vec2(20.0, 14.0)), 3.0, chip);
            if *parcel {
                painter.circle(p, dot_r.clamp(3.0, 4.5), *c, Stroke::new(0.8, Color32::from_black_alpha(140)));
            } else {
                // Spark with a short zig-zag trail, as drawn on the faces.
                let trail = [p, p + vec2(-3.0, 2.5), p + vec2(-6.0, -1.0), p + vec2(-9.0, 2.0)];
                for w in trail.windows(2) {
                    painter.line_segment([w[0], w[1]], Stroke::new(1.0, PHOTON.gamma_multiply(0.5)));
                }
                photon(painter, p, 3.0, 1.0);
            }
        }
        painter.text(
            pos2(tx, y),
            Align2::LEFT_BOTTOM,
            text,
            font.clone(),
            if icon.is_some() { ink } else { secondary },
        );
        y -= 18.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn octant_patch_is_a_closed_triangulation() {
        let (mut verts, mut tris) = (Vec::new(), Vec::new());
        octant_patch(&mut verts, &mut tris, [1.0, -1.0, 1.0], 0.5, [1.0; 3], Kind::Core);
        assert_eq!(verts.len(), (N + 1) * (N + 2) / 2);
        assert_eq!(tris.len(), N * N);
        assert!(tris.iter().flatten().all(|&i| (i as usize) < verts.len()));
        assert!(verts.iter().all(|v| (dot(v.p, v.p).sqrt() - 0.5).abs() < 1e-5));
    }

    #[test]
    fn visibility_of_cut_faces() {
        let view_into_cut = normalize([1.0, 1.0, 1.0]);
        // A point in the middle of face z = 0, seen from inside the open octant.
        let p = face_point(0, 0.6, 0.7);
        assert!(visible(p, view_into_cut, 0.0));
        // The same point seen from below is hidden by the solid material.
        assert!(!visible(p, [0.0, 0.0, -1.0], 0.0));
        // Seen from a grazing direction through another octant: hidden.
        assert!(!visible(p, normalize([-1.0, 0.1, 0.3]), 0.0));
        // A face point (outside the core ball) whose line of sight hits the
        // ball's visible octant at (0.2, 0.1, 0.15) is hidden by the ball.
        let p = [0.35, 0.05, 0.0];
        let d = normalize([-0.15, 0.05, 0.15]);
        assert!(!visible(p, d, 0.3));
        // Looking straight out of the face, past the ball: visible.
        assert!(visible(p, [0.0, 0.0, 1.0], 0.3));
    }

    #[test]
    fn ring_radii_increase_and_end_at_surface() {
        let rn: Vec<f64> = (0..=500).map(|i| 1.0 - (1.0 - i as f64 / 500.0).powi(4)).collect();
        let r = ring_radii(&rn, 0.2, &[0.7]);
        assert_eq!(r[0], 0.2);
        assert_eq!(*r.last().unwrap(), 1.0);
        assert!(r.windows(2).all(|w| w[1] > w[0]));
        assert!(r.iter().any(|&x| (x - 0.7).abs() < 1e-6));
    }
}
