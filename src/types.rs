#[derive(Clone, Copy, Debug)]
pub struct ComplexDouble {
    pub re: f64,
    pub im: f64,
}

impl ComplexDouble {
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub fn norm(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }

    pub fn sqrt(self) -> Self {
        let magnitude = self.norm();
        if magnitude == 0.0 {
            return Self::new(0.0, 0.0);
        }
        let real_part = ((magnitude + self.re) / 2.0).sqrt();
        let imag_part = if self.im >= 0.0 {
            self.im.abs() / (2.0 * real_part)
        } else {
            -(self.im.abs() / (2.0 * real_part))
        };
        Self::new(real_part, imag_part)
    }

    pub fn exp(self) -> Self {
        let exp_re = self.re.exp();
        Self::new(exp_re * self.im.cos(), exp_re * self.im.sin())
    }
}

impl std::ops::Add for ComplexDouble {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl std::ops::Sub for ComplexDouble {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl std::ops::Mul for ComplexDouble {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self::new(self.re * rhs.re - self.im * rhs.im, self.re * rhs.im + self.im * rhs.re)
    }
}

impl std::ops::Div for ComplexDouble {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let denom = rhs.re * rhs.re + rhs.im * rhs.im;
        Self::new((self.re * rhs.re + self.im * rhs.im) / denom, (self.im * rhs.re - self.re * rhs.im) / denom)
    }
}

impl std::ops::Mul<f64> for ComplexDouble {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.re * rhs, self.im * rhs)
    }
}

impl std::ops::Neg for ComplexDouble {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}

#[derive(Clone, Debug)]
pub struct IntermediateValues {
    pub theta_hzn: [f64; 2],
    pub d_hzn__meter: [f64; 2],
    pub h_e__meter: [f64; 2],
    pub n_s: f64,
    pub delta_h__meter: f64,
    pub a_ref__db: f64,
    pub a_fs__db: f64,
    pub d__km: f64,
    pub mode: i32,
}

impl Default for IntermediateValues {
    fn default() -> Self {
        Self {
            theta_hzn: [0.0; 2],
            d_hzn__meter: [0.0; 2],
            h_e__meter: [0.0; 2],
            n_s: 0.0,
            delta_h__meter: 0.0,
            a_ref__db: 0.0,
            a_fs__db: 0.0,
            d__km: 0.0,
            mode: 0,
        }
    }
}