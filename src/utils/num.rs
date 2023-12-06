

pub fn range_scale(v: f32, old_low: f32, old_hi: f32, new_low: f32, new_hi: f32) -> f32 {
    // Scale 'v' from ['old_low', 'old_hi'] to ['new_low', 'new_hi']
    return new_low + v * (new_hi - new_low) / (old_hi - old_low);
}

pub fn clamp(val: f32, min: f32, max: f32) -> f32 {
    // Clamp 'val' between 'min' and 'max'
    match val {
        val if val < min => min,
        val if val > max => max,
        _ => val
    }
}

pub fn rvec2_range(rng: &mut ThreadRng, x_min: f32, x_max: f32, y_min: f32, y_max: f32) -> Vector2 {
    rvec2(range_scale(rng.gen::<f32>(), 0., 1., x_min, x_max),
        range_scale(rng.gen::<f32>(), 0., 1., y_min, y_max))
}