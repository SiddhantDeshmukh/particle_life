use std::f32::consts::PI;

use rand::{rngs::ThreadRng, Rng};
use raylib::prelude::*;
use rayon::prelude::*;

// Defined for convenience, needs refactor
pub const X_SIZE: f32 = 100.;
pub const Y_SIZE: f32 = 100.;

pub struct Button {
    pub x: f32,
    pub y: f32,
    pub width: f32,
    pub height: f32,
    pub text: str,
}


#[derive(Debug)]
pub enum BoundaryCondition {
    Periodic,
    Reflecting
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub struct Particle {
    pub color: Color,
    pub i_color: usize,  // index of Color in colors
    pub color_vec: [f32; 3],  // Color saved as RGB vec for continuous force
    pub position: Vector2,
    pub velocity: Vector2,
}

impl Particle {
    // Scale positions to window size for drawing, divide by domain size
    fn scale_position(&self, value: f32, window_size: i32) -> i32 {
        (value * window_size as f32) as i32
    }

    pub fn win_pos_x(&self, window_width: i32) -> i32 {
        self.scale_position(self.position.x, window_width) / X_SIZE as i32
    }

    pub fn win_pos_y(&self, window_height: i32) -> i32 {
        self.scale_position(self.position.y, window_height) / Y_SIZE as i32
    }
}

// struct SpawnBounds {
//     // Controls particle spawning
//     x_min: f32,
//     x_max: f32,
//     y_min: f32,
//     y_max: f32,
//     vx_min: f32,
//     vx_max: f32,
//     vy_min: f32,
//     vy_max: f32
// }

// impl SpawnBounds {
//     fn full_vzero(params: Params) -> Self {
//         // Spawn bounds is the entire window, velocities are zero
//         Self {
//             x_min: 0.,
//             x_max: params.width as f32,
//             y_min: 0.,
//             y_max: params.height as f32,
//             vx_min: 0.,
//             vx_max: 0.,
//             vy_min: 0.,
//             vy_max: 0.,
//         }
//     }
// }


#[derive(Debug)]
pub struct Params {
    pub window_width: i32,
    pub window_height: i32,
    pub time_step: f32,
    pub friction_half_life: f32,
    pub max_radius: f32,
    pub boundary_condition: BoundaryCondition,
    pub grid_prgb_arr: [f32; 3],  // 1-way force from grid cell to particle
    pub force_scale: f32,
}

impl Params {
    pub fn x_min(&self) -> f32 {
        0.  // default
    }

    pub fn x_max(&self) -> f32 {
        X_SIZE  // default        
    }

    pub fn y_min(&self) -> f32 {
        0.  // default
    }

    pub fn y_max(&self) -> f32 {
        Y_SIZE  // default
    }

    pub fn x_len(&self) -> f32 {
        self.x_max() - self.x_min()
    }

    pub fn y_len(&self) -> f32 {
        self.y_max() - self.y_min()
    }
}

// Coordinate transforms
pub fn win_to_world(vec: Vector2, window_width: i32, window_height: i32) -> Vector2 {
    // Convert from window coordinates to world coords (between 0 and 1)
    Vector2{x: vec.x * X_SIZE / window_width as f32, y: vec.y * Y_SIZE / window_height as f32}
}

pub fn world_to_win(vec: Vector2, window_width: i32, window_height: i32) -> Vector2 {
    Vector2{x: vec.x / X_SIZE * window_width as f32, y: vec.y / Y_SIZE * window_height as f32}
}

// Numerics
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

pub fn random_val(rng: &mut ThreadRng) -> f32 {
    // Convenience method
    (rng.gen::<f32>() - 0.5) * 2.
}

// Force matrix stuff

pub fn generate_prgb_matrix(rng: &mut ThreadRng) -> [f32; 3] {
    // 3 random values in range [-1., 1.] for continuous color force
    [random_val(rng), random_val(rng), random_val(rng)]
}

pub fn generate_force_matrix(len: usize, rng: &mut ThreadRng) -> Vec<Vec<f32>> {
    // Generate square matrix init with random values in [-1., 1.]
    let mut force_matrix = vec![vec![0 as f32; len]; len];
    // Populate
    for i in 0..len {
        for j in 0..len {
            force_matrix[i][j] = (rng.gen::<f32>() - 0.5) * 2.;
        }
    }

    return force_matrix
}

pub fn compute_force(r: f32, a: f32, beta: f32) -> f32 {
    match r {
        r if r < beta => r / beta - 1.0,
        r if beta < r && r < 1.0 => a * (1.0 - (2.0 * r - 1.0 - beta).abs()) / (1.0 - beta),
        _ => 0.0,
    }
}

// Color forces
pub fn generate_rgb_matrix(rng: &mut ThreadRng) -> Vec<Vec<f32>> {
    // A random 3x3 matrix represent (RGB) cross-correlations
    generate_force_matrix(3, rng)
}

pub fn random_color_vec(rng: &mut ThreadRng) -> [f32; 3] {
    let mut color_vec = [0., 0., 0.];
    for i in 0..3 {
        color_vec[i] = rng.gen_range(0..255) as f32;
    }
    color_vec
}

pub fn color_attraction(cv1: [f32; 3], cv2: [f32; 3], rgb_matrix: &Vec<Vec<f32>>) -> f32 {
    // Continuous force calculation between two colors, designed to be
    // asymmetric, but similar colors exert smaller forces
    let mut force: f32 = 0.;
    for i in 0..3 {
        for j in 0..3 {
            // sin of L1 norm
            let a: f32 = range_scale(cv1[i] - cv2[j], -255., 255., 0., 2. * PI).sin();
            // Random number in range [-1., 1.]
            let b = rgb_matrix[i][j];
            force += a + b;
        }
    }
    return force.cos();
}

pub fn color_to_vec(c: Color) -> [f32; 3] {
    // Turn a Color into a (RGB) vec, ignoring the alppha component
    [c.r as f32, c.g as f32, c.b as f32]
}

pub fn vec_to_color(v: [f32; 3]) -> Color {
    // Turn an RGB vec into a Color, alpha always max
    Color{r: v[0] as u8, g: v[1] as u8, b: v[2] as u8, a: 255}
}

pub fn random_color(rng: &mut ThreadRng) -> Color {
    vec_to_color(random_color_vec(rng))
}

// Particle generation
pub fn generate_particles_cm(num: usize, rng: &mut ThreadRng,
                     colors: &Vec<Color>, x_min: f32, x_max: f32,
                     y_min: f32, y_max: f32, vx_min: f32, vx_max: f32,
                     vy_min: f32, vy_max: f32) -> Vec<Particle> {
    // Generate 'num' random particles from the given parameters
    let particles: Vec<Particle> = (0..num)
        .map(|_|  {
            let i_color: usize = rng.gen_range(0..colors.len());
            let p = Particle {
                color: colors[i_color],
                i_color,
                color_vec: color_to_vec(colors[i_color]),  // used for grid forces
                position: rvec2_range(rng, x_min, x_max, y_min, y_max),
                velocity: rvec2_range(rng, vx_min, vx_max, vy_min, vy_max),
            };
            p
        })
        .collect();
    particles
}

// Particle updates
pub fn update_particles_cm(particles: &Vec<Particle>, params: &Params,
                    force_matrix: &Vec<Vec<f32>>) -> Vec<Particle> {
    // Update using a fixed color matrix
    let friction: f32 = 0.5_f32.powf(params.time_step / params.friction_half_life);
    // Create new particles list and update velocities
    let new_particles: Vec<Particle> = particles
        .par_iter()
        .map(|p1| {
            let mut total_force = Vector2::zero();
            for p2 in particles {
                if p1 == p2 {
                    // No self-attraction
                    continue;
                };
                let distance: f32;
                let mut p2_pos = p2.position;
                if let BoundaryCondition::Periodic = params.boundary_condition {
                    // Find closest periodic image of the target point
                    let distance_x = p2_pos.x - p1.position.x;
                    if distance_x.abs() > 0.5 * params.x_len() {
                        if distance_x > 0. {
                            // Nearest image is to the left
                            p2_pos.x -= params.x_len();
                        } else {
                            // Nearest image is to the right
                            p2_pos.x += params.x_len();
                        }
                    }

                    let distance_y = p2_pos.y - p1.position.y;
                    if distance_y.abs() > 0.5 * params.y_len() {
                        if distance_y > 0. {
                            // Nearest image is below
                            p2_pos.y -= params.y_len();
                        } else {
                            // Nearest image is above
                            p2_pos.y += params.y_len();
                        }
                    }
                };
                distance = p1.position.distance_to(p2_pos);
                if (distance > 0.) & (distance < params.max_radius) {
                    let f = compute_force(
                        distance / params.max_radius,
                        force_matrix[p1.i_color][p2.i_color],
                        0.3,
                    );
                    total_force += ((p2_pos - p1.position) / distance) * f;
                }
            }
            let mut new_p = *p1;
            new_p.velocity *= friction;
            total_force *= params.force_scale / (X_SIZE * Y_SIZE);
            new_p.velocity += total_force * params.time_step;
            new_p
        })
        .collect();

    // Update positions
    new_particles
        .par_iter()
        .map(|p| {
            let mut new_p = *p;
            // Potential new position
            new_p.position += new_p.velocity * params.time_step;
            // Boundaries
            match params.boundary_condition {
                BoundaryCondition::Periodic => {
                    match new_p.position.x {
                        x if x > params.x_max() => {
                            new_p.position.x -= params.x_len();
                        },
                        x if x < params.x_min() => {
                            new_p.position.x += params.x_len();
                        },
                        _ => {}
                    };
                    match new_p.position.y {
                        y if y > params.y_max() => {
                            new_p.position.y -= params.y_len();
                        },
                        y if y < params.y_min() => {
                            new_p.position.y += params.y_len();
                        },
                        _ => {}
                    }
                },
                BoundaryCondition::Reflecting => {
                    match new_p.position.x {
                        x if x > params.x_max() =>  {
                            new_p.position.x = 2. * params.x_max() - new_p.position.x;
                            new_p.velocity.x = -new_p.velocity.x;
                        },
                        x if x < params.x_min() => {
                            new_p.position.x = 2. * params.x_min() - new_p.position.x;
                            new_p.velocity.x = -new_p.velocity.x;
                        },
                        _ => {}
                    };
                    match new_p.position.y {
                        y if y > params.y_max() => {
                            new_p.position.y = 2. * params.y_max() - new_p.position.y;
                            new_p.velocity.y = -new_p.velocity.y;
                        },
                        y if y < params.y_min() => {
                            new_p.position.y = 2. * params.y_min() - new_p.position.y;
                            new_p.velocity.y = -new_p.velocity.y;
                        },
                        _ => {}
                    };
               },
            }
            new_p
        })
        .collect()
}