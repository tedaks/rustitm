pub fn itm_max(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}

pub fn itm_min(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}

pub fn fortran_dim(x: f64, y: f64) -> f64 {
    if x > y {
        x - y
    } else {
        0.0
    }
}

pub fn fresnel_integral(v2: f64) -> f64 {
    if v2 < 5.76 {
        6.02 + 9.11 * v2.sqrt() - 1.27 * v2
    } else {
        12.953 + 10.0 * v2.log10()
    }
}

pub fn height_function(x__km: f64, k: f64) -> f64 {
    let w: f64;
    let mut result: f64;

    if x__km < 200.0 {
        w = -k.ln();

        if k < 1e-5 || x__km * w.powi(3) > 5495.0 {
            result = -117.0;
            if x__km > 1.0 {
                result = 17.372 * x__km.ln() + result;
            }
        } else {
            result = 2.5e-5 * x__km.powi(2) / k - 8.686 * w - 15.0;
        }
    } else {
        result = 0.05751 * x__km - 4.343 * x__km.ln();

        if x__km < 2000.0 {
            w = 0.0134 * x__km * (-0.005 * x__km).exp();
            result = (1.0 - w) * result + w * (17.372 * x__km.ln() - 117.0);
        }
    }

    result
}

pub fn h0_curve(j: i32, r: f64) -> f64 {
    let a = [25.0, 80.0, 177.0, 395.0, 705.0];
    let b = [24.0, 45.0, 68.0, 80.0, 105.0];

    10.0 * (1.0 + a[j as usize] * (1.0 / r).powi(4) + b[j as usize] * (1.0 / r).powi(2)).log10()
}

pub fn h0_function(r: f64, mut eta_s: f64) -> f64 {
    eta_s = eta_s.max(1.0).min(5.0);

    let i = eta_s as i32;
    let q = eta_s - i as f64;

    let mut result = h0_curve(i - 1, r);

    if q != 0.0 {
        result = (1.0 - q) * result + q * h0_curve(i, r);
    }

    result
}

pub fn inverse_complementary_cumulative_distribution_function(q: f64) -> f64 {
    let c_0 = 2.515516;
    let c_1 = 0.802853;
    let c_2 = 0.010328;
    let d_1 = 1.432788;
    let d_2 = 0.189269;
    let d_3 = 0.001308;

    let mut x = q;
    if q > 0.5 {
        x = 1.0 - x;
    }

    let t_x = (-2.0 * x.ln()).sqrt();
    let zeta_x = ((c_2 * t_x + c_1) * t_x + c_0) / (((d_3 * t_x + d_2) * t_x + d_1) * t_x + 1.0);

    let mut q_q = t_x - zeta_x;

    if q > 0.5 {
        q_q = -q_q;
    }

    q_q
}

pub fn curve(c1: f64, c2: f64, x1: f64, x2: f64, x3: f64, d_e__meter: f64) -> f64 {
    (c1 + c2 / (1.0 + ((d_e__meter - x2) / x3).powi(2))) * (d_e__meter / x1).powi(2)
        / (1.0 + (d_e__meter / x1).powi(2))
}
