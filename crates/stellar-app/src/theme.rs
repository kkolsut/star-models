//! Colors: the validated categorical palette (light/dark steps), neutral chart
//! ink, and single-hue sequential ramps for the cross-section.
//!
//! Categorical order is fixed and never cycled. The current model always uses
//! slot 0; each pinned model keeps the slot it received when it was pinned.

use egui::Color32;

const fn hex(v: u32) -> Color32 {
    Color32::from_rgb((v >> 16) as u8, (v >> 8) as u8, v as u8)
}

/// Reference categorical palette (validated: adjacent CVD ΔE ≥ 8, normal-vision ΔE ≥ 15).
pub const SERIES_LIGHT: [Color32; 8] = [
    hex(0x2a78d6), // blue
    hex(0xeb6834), // orange
    hex(0x1baf7a), // aqua
    hex(0xeda100), // yellow
    hex(0xe87ba4), // magenta
    hex(0x008300), // green
    hex(0x6250d6), // violet
    hex(0xe34948), // red
];
pub const SERIES_DARK: [Color32; 8] = [
    hex(0x3987e5),
    hex(0xd95926),
    hex(0x199e70),
    hex(0xc98500),
    hex(0xd55181),
    hex(0x008300),
    hex(0x9085e9),
    hex(0xe66767),
];
pub const MAX_SLOTS: usize = 8;

pub fn series(slot: usize, dark: bool) -> Color32 {
    let pal = if dark { &SERIES_DARK } else { &SERIES_LIGHT };
    pal[slot.min(MAX_SLOTS - 1)]
}

/// Muted ink for axes, reference lines (ZAMS) and annotations.
pub fn muted(_dark: bool) -> Color32 {
    hex(0x898781)
}

/// Secondary ink for text that should recede a little.
pub fn secondary(dark: bool) -> Color32 {
    if dark { hex(0xc3c2b7) } else { hex(0x52514e) }
}

/// Primary ink.
pub fn primary(dark: bool) -> Color32 {
    if dark { hex(0xf0efec) } else { hex(0x0b0b0b) }
}

/// Neutral translucent fill marking convective zones (never a series color).
pub fn convective_fill(dark: bool) -> Color32 {
    if dark {
        Color32::from_rgba_unmultiplied(195, 194, 183, 26)
    } else {
        Color32::from_rgba_unmultiplied(137, 135, 129, 40)
    }
}

// ---------------------------------------------------------------------------
// OKLab helpers and sequential ramps

fn srgb_to_linear(c: u8) -> f64 {
    let c = c as f64 / 255.0;
    if c <= 0.04045 { c / 12.92 } else { ((c + 0.055) / 1.055).powf(2.4) }
}

fn linear_to_srgb(c: f64) -> f64 {
    if c <= 0.0031308 { 12.92 * c } else { 1.055 * c.powf(1.0 / 2.4) - 0.055 }
}

fn to_oklab(c: Color32) -> [f64; 3] {
    let (r, g, b) = (srgb_to_linear(c.r()), srgb_to_linear(c.g()), srgb_to_linear(c.b()));
    let l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b).cbrt();
    let m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b).cbrt();
    let s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b).cbrt();
    [
        0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
        1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
        0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s,
    ]
}

/// OKLab → linear sRGB (may be out of gamut).
fn oklab_to_linear(lab: [f64; 3]) -> [f64; 3] {
    let [l, a, b] = lab;
    let l_ = (l + 0.3963377774 * a + 0.2158037573 * b).powi(3);
    let m_ = (l - 0.1055613458 * a - 0.0638541728 * b).powi(3);
    let s_ = (l - 0.0894841775 * a - 1.2914855480 * b).powi(3);
    [
        4.0767416621 * l_ - 3.3077115913 * m_ + 0.2309699292 * s_,
        -1.2684380046 * l_ + 2.6097574011 * m_ - 0.3413193965 * s_,
        -0.0041960863 * l_ - 0.7034186147 * m_ + 1.7076147010 * s_,
    ]
}

fn from_linear(rgb: [f64; 3]) -> Color32 {
    let q = |v: f64| (linear_to_srgb(v.clamp(0.0, 1.0)) * 255.0).round() as u8;
    Color32::from_rgb(q(rgb[0]), q(rgb[1]), q(rgb[2]))
}

/// OKLCH color, with chroma reduced until it fits the sRGB gamut.
fn oklch(l: f64, c: f64, h: f64) -> Color32 {
    let mut c = c;
    loop {
        let rgb = oklab_to_linear([l, c * h.cos(), c * h.sin()]);
        if rgb.iter().all(|v| (-1e-6..=1.0 + 1e-6).contains(v)) || c < 1e-4 {
            return from_linear(rgb);
        }
        c *= 0.97;
    }
}

/// A 256-entry color lookup table, ordered from low to high magnitude for the
/// active theme. The low end recedes toward the chart surface: light in light
/// mode, dark in dark mode.
#[derive(Clone)]
pub struct Ramp {
    lut: Vec<Color32>,
}

impl Ramp {
    pub fn sample(&self, t: f64) -> Color32 {
        let t = if t.is_finite() { t.clamp(0.0, 1.0) } else { 0.0 };
        self.lut[(t * (self.lut.len() - 1) as f64).round() as usize]
    }

    fn oriented(mut lut: Vec<Color32>, dark: bool) -> Self {
        // `lut` runs light → dark; dark mode puts low values at the dark end.
        if dark {
            lut.reverse();
        }
        Self { lut }
    }

    /// Reference blue sequential ramp (steps 100–700), interpolated in OKLab.
    pub fn blue(dark: bool) -> Self {
        const STEPS: [u32; 13] = [
            0xcde2fb, 0xb7d3f6, 0x9ec5f4, 0x86b6ef, 0x6da7ec, 0x5598e7, 0x3987e5, 0x2a78d6, 0x256abf, 0x1c5cab,
            0x184f95, 0x104281, 0x0d366b,
        ];
        let labs: Vec<[f64; 3]> = STEPS.iter().map(|&v| to_oklab(hex(v))).collect();
        let lut = (0..256)
            .map(|i| {
                let x = i as f64 / 255.0 * (labs.len() - 1) as f64;
                let k = (x.floor() as usize).min(labs.len() - 2);
                let w = x - k as f64;
                let mix = |j: usize| labs[k][j] + w * (labs[k + 1][j] - labs[k][j]);
                from_linear(oklab_to_linear([mix(0), mix(1), mix(2)]))
            })
            .collect();
        Self::oriented(lut, dark)
    }

    /// Single-hue orange ramp at the hue of categorical slot 2, from
    /// L = 0.93 down to L = 0.36 in OKLCH.
    pub fn orange(dark: bool) -> Self {
        let [_, a, b] = to_oklab(hex(0xeb6834));
        let h = b.atan2(a);
        let c_peak = (a * a + b * b).sqrt();
        let lut = (0..256)
            .map(|i| {
                let t = i as f64 / 255.0;
                let l = 0.93 - 0.57 * t;
                let c = c_peak * (1.0 - ((l - 0.62) / 0.42).powi(2)).max(0.25);
                oklch(l, c, h)
            })
            .collect();
        Self::oriented(lut, dark)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lightness(c: Color32) -> f64 {
        to_oklab(c)[0]
    }

    #[test]
    fn ramps_are_monotonic_in_lightness() {
        for ramp in [Ramp::blue(false), Ramp::orange(false)] {
            let ls: Vec<f64> = ramp.lut.iter().map(|&c| lightness(c)).collect();
            assert!(ls.windows(2).all(|w| w[1] <= w[0] + 2e-3), "light mode must darken with magnitude");
        }
        for ramp in [Ramp::blue(true), Ramp::orange(true)] {
            let ls: Vec<f64> = ramp.lut.iter().map(|&c| lightness(c)).collect();
            assert!(ls.windows(2).all(|w| w[1] >= w[0] - 2e-3), "dark mode must lighten with magnitude");
        }
    }

    #[test]
    fn oklab_round_trip() {
        for &v in &[0x2a78d6, 0xeb6834, 0x0b0b0b, 0xfcfcfb] {
            let c = hex(v);
            assert_eq!(from_linear(oklab_to_linear(to_oklab(c))), c);
        }
    }
}
