use crate::constants::PI;
use crate::helper::{itm_max, itm_min};
use crate::terrain::{sigma_h_function, terrain_roughness};
use crate::types::ComplexDouble;

pub fn line_of_sight_loss(
    d__meter: f64,
    h_e__meter: [f64; 2],
    z_g: ComplexDouble,
    delta_h__meter: f64,
    m_d: f64,
    a_d0: f64,
    d_sML__meter: f64,
    f__mhz: f64,
) -> f64 {
    let delta_h_d__meter = terrain_roughness(d__meter, delta_h__meter);
    let sigma_h_d__meter = sigma_h_function(delta_h_d__meter);

    let wn = f__mhz / 47.7;

    let sin_psi = (h_e__meter[0] + h_e__meter[1])
        / (d__meter.powi(2) + (h_e__meter[0] + h_e__meter[1]).powi(2)).sqrt();

    let sin_psi_c = ComplexDouble::new(sin_psi, 0.0);
    let mut r_e = (sin_psi_c - z_g) / (sin_psi_c + z_g)
        * (-itm_min(10.0, wn * sigma_h_d__meter * sin_psi)).exp();

    let q = r_e.re * r_e.re + r_e.im * r_e.im;
    if q < 0.25 || q < sin_psi {
        r_e = r_e * (sin_psi / q).sqrt();
    }

    let mut delta_phi = wn * 2.0 * h_e__meter[0] * h_e__meter[1] / d__meter;

    if delta_phi > PI / 2.0 {
        delta_phi = PI - (PI / 2.0).powi(2) / delta_phi;
    }

    let rr = ComplexDouble::new(delta_phi.cos(), -delta_phi.sin()) + r_e;
    let a_t__db = -10.0 * (rr.re.powi(2) + rr.im.powi(2)).log10();

    let a_d__db = m_d * d__meter + a_d0;

    let w = 1.0 / (1.0 + f__mhz * delta_h__meter / itm_max(10e3, d_sML__meter));

    w * a_t__db + (1.0 - w) * a_d__db
}
