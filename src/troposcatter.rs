use crate::constants::SQRT2;
use crate::helper::{h0_function, itm_max, itm_min};

pub fn f_function(td: f64) -> f64 {
    let a = [133.4, 104.6, 71.8];
    let b = [0.332e-3, 0.212e-3, 0.157e-3];
    let c = [-10.0, -2.5, 5.0];

    let i = if td <= 10e3 {
        0
    } else if td <= 70e3 {
        1
    } else {
        2
    };

    a[i] + b[i] * td + c[i] * td.log10()
}

pub fn troposcatter_loss(
    d__meter: f64,
    theta_hzn: [f64; 2],
    d_hzn__meter: [f64; 2],
    h_e__meter: [f64; 2],
    a_e__meter: f64,
    n_s: f64,
    f__mhz: f64,
    theta_los: f64,
    h0: &mut f64,
) -> f64 {
    let wn = f__mhz / 47.7;

    let h_0_val;

    if *h0 > 15.0 {
        h_0_val = *h0;
    } else {
        let mut ad = d_hzn__meter[0] - d_hzn__meter[1];
        let mut rr = h_e__meter[1] / h_e__meter[0];

        if ad < 0.0 {
            ad = -ad;
            rr = 1.0 / rr;
        }

        let theta = theta_hzn[0] + theta_hzn[1] + d__meter / a_e__meter;

        let r_1 = 2.0 * wn * theta * h_e__meter[0];
        let r_2 = 2.0 * wn * theta * h_e__meter[1];

        if r_1 < 0.2 && r_2 < 0.2 {
            return 1001.0;
        }

        let mut s = (d__meter - ad) / (d__meter + ad);
        let q = itm_min(itm_max(0.1, rr / s), 10.0);
        s = itm_max(0.1, s);

        let h_0__meter = (d__meter - ad) * (d__meter + ad) * theta * 0.25 / d__meter;

        let z_0__meter = 1.7556e3;
        let z_1__meter = 8.0e3;
        let eta_s = (h_0__meter / z_0__meter)
            * (1.0
                + (0.031 - n_s * 2.32e-3 + n_s.powi(2) * 5.67e-6)
                    * (-itm_min(1.7, h_0__meter / z_1__meter).powi(6)).exp());

        let h_00 = (h0_function(r_1, eta_s) + h0_function(r_2, eta_s)) / 2.0;
        let delta_h_0 = itm_min(
            h_00,
            6.0 * (0.6 - itm_max(eta_s, 1.0).log10()) * s.log10() * q.log10(),
        );

        let mut h_0_result = h_00 + delta_h_0;
        h_0_result = itm_max(h_0_result, 0.0);

        if eta_s < 1.0 {
            h_0_result = eta_s * h_0_result
                + (1.0 - eta_s)
                    * 10.0
                    * (((1.0 + SQRT2 / r_1) * (1.0 + SQRT2 / r_2)).powi(2) * (r_1 + r_2)
                        / (r_1 + r_2 + 2.0 * SQRT2))
                        .log10();
        }

        if h_0_result > 15.0 && *h0 >= 0.0 {
            h_0_result = *h0;
        }

        h_0_val = h_0_result;
    }

    *h0 = h_0_val;
    let th = d__meter / a_e__meter - theta_los;

    let d_0__meter = 40e3;
    let h__meter = 47.7;
    f_function(th * d__meter) + 10.0 * (wn * h__meter * th.powi(4)).log10()
        - 0.1 * (n_s - 301.0) * (-th * d__meter / d_0__meter).exp()
        + h_0_val
}
