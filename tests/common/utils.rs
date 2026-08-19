// SPDX-License-Identifier: LGPL-3.0-or-later

#[allow(dead_code)]
pub fn modinv_i128(a: i128, m: i128) -> Option<i128> {
    let (mut t, mut new_t) = (0, 1);
    let (mut r, mut new_r) = (m, a);
    while new_r != 0 {
        let q = r / new_r;
        (t, new_t) = (new_t, t - q * new_t);
        (r, new_r) = (new_r, r - q * new_r);
    }
    if r > 1 {
        return None;
    }
    if t < 0 {
        t += m;
    }
    return Some(t);
}
