use rand::{rngs::ThreadRng, Rng};
use raylib::prelude::*;
use rayon::prelude::*;

enum BoundaryCondition {
    Periodic,
    Reflecting
}

#[derive(PartialEq, Clone, Copy)]
struct Particle {
    color: Color,
    i_color: usize,  // index of Color in colors
    position: Vector2,
    velocity: Vector2,
    mass: f32
}

struct SpawnBounds {
    // Controls particle spawning
    x_min: f32,
    x_max: f32,
    y_min: f32,
    y_max: f32,
    vx_min: f32,
    vx_max: f32,
    vy_min: f32,
    vy_max: f32
}

impl SpawnBounds {
    fn full_vzero(params: Params) -> Self {
        // Spawn bounds is the entire window, velocities are zero
        Self {
            x_min: 0.,
            x_max: params.width as f32,
            y_min: 0.,
            y_max: params.height as f32,
            vx_min: 0.,
            vx_max: 0.,
            vy_min: 0.,
            vy_max: 0.,
        }
    }
}

struct Params {
    // width and height are for raylib
    width: i32,
    height: i32,
    friction_half_life: f32,
    time_step: f32,
    max_radius: f32,
    boundary_condition: BoundaryCondition,
}

impl Params {
    fn x_min(&self) -> f32 {
        0.  // default
    }

    fn x_max(&self) -> f32 {
        self.width as f32
    }

    fn y_min(&self) -> f32 {
        0.  // default
    }

    fn y_max(&self) -> f32 {
        self.height as f32
    }

    fn x_len(&self) -> f32 {
        self.x_max() - self.x_min()
    }

    fn y_len(&self) -> f32 {
        self.y_max() - self.y_min()
    }
}

fn generate_force_matrix(len: usize, rng: &mut ThreadRng) -> Vec<Vec<f32>> {
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

fn range_scale(v: f32, old_low: f32, old_hi: f32, new_low: f32, new_hi: f32) -> f32 {
    // Scale 'v' from ['old_low', 'old_hi'] to ['new_low', 'new_hi']
    return new_low + v * (new_hi - new_low) / (old_hi - old_low);
}

fn rvec2_range(rng: &mut ThreadRng, x_min: f32, x_max: f32, y_min: f32, y_max: f32) -> Vector2 {
    rvec2(range_scale(rng.gen::<f32>(), 0., 1., x_min, x_max),
          range_scale(rng.gen::<f32>(), 0., 1., y_min, y_max))
}

fn compute_force(r: f32, a: f32, beta: f32) -> f32 {
    match r {
        r if r < beta => r / beta - 1.0,
        r if beta < r && r < 1.0 => a * (1.0 - (2.0 * r - 1.0 - beta).abs()) / (1.0 - beta),
        _ => 0.0,
    }
}

fn generate_particles(num: usize, rng: &mut ThreadRng,
                     colors: Vec<Color>, x_min: f32, x_max: f32,
                     y_min: f32, y_max: f32, vx_min: f32, vx_max: f32,
                     vy_min: f32, vy_max: f32) -> Vec<Particle> {
    // Generate 'num' random particles from the given parameters
    let particles: Vec<Particle> = (0..num)
        .map(|_|  {
            let i_color: usize = rng.gen_range(0..colors.len());
            let p = Particle {
                color: colors[i_color],
                i_color,
                position: rvec2_range(rng, x_min, x_max, y_min, y_max),
                velocity: rvec2_range(rng, vx_min, vx_max, vy_min, vy_max),
                mass: rng.gen_range(1.0..8.0),
            };

            p
        })
        .collect();
    return particles;
}

fn update_particles(particles: &Vec<Particle>, params: &Params,
                    force_matrix: &Vec<Vec<f32>>) -> Vec<Particle> {
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
                if let BoundaryCondition::Periodic = params.boundary_condition {
                    // Check component-wise if distances exceed the half-domain and adjust if necessary
                    let mut pos_change = p2.position - p1.position;
                    pos_change.x = pos_change.x.abs();
                    pos_change.y = pos_change.y.abs();
                    if pos_change.x > 0.5 * params.x_len() {
                        pos_change.x = pos_change.x - params.x_len();
                    }
                    if pos_change.y > 0.5 * params.y_len() {
                        pos_change.y = pos_change.y - params.y_len();
                    }
                    distance = pos_change.length();
                    // distance = p1.position.distance_to(p2.position);
                } else {
                    distance = p1.position.distance_to(p2.position);
                };
                if (distance > 0.) & (distance < params.max_radius) {
                    let f = compute_force(
                        distance / params.max_radius,
                        force_matrix[p1.i_color][p2.i_color],
                        0.3,
                    );
                    total_force += ((p2.position - p1.position) / distance) * f * p2.mass;
                }
            }
            let mut new_p = *p1;
            new_p.velocity *= friction;
            new_p.velocity += total_force * params.time_step;

            // Boundaries
            match params.boundary_condition {
                BoundaryCondition::Periodic => {
                    match new_p.position.x {
                        x if x > params.x_max() => {
                            new_p.position.x -= params.x_max();
                        },
                        x if x < params.x_min() => {
                            new_p.position.x += params.x_max();
                        },
                        _ => {}
                    };
                    match new_p.position.y {
                        y if y > params.y_max() => {
                            new_p.position.y -= params.y_max();
                        },
                        y if y < params.y_min() => {
                            new_p.position.y += params.y_max();
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
        .collect();

    // Update positions
    new_particles
        .iter()
        .map(|p| {
            let mut new_p = *p;
            new_p.position += new_p.velocity * params.time_step;
            new_p
        })
        .collect()
}

fn main() {
    // Simulation parameters
    let params = Params {
        width: 1200,
        height: 800,
        friction_half_life: 0.1,
        time_step: 0.05,
        max_radius: 200.,
        boundary_condition: BoundaryCondition::Periodic
        // boundary_condition: BoundaryCondition::Reflecting
    };

    // let window_padding = 100;

    const NUM_PARTICLES: usize = 1000;
    // const SEED: u64 = 420;

    // Bounds for initial particle spawning
    let x_min: f32 = 0.;
    let x_max: f32 = params.width as f32;
    let y_min: f32 = 0.;
    let y_max: f32 = params.height as f32;
    let vx_min: f32 = 0.;
    let vx_max: f32 = 0.;
    let vy_min: f32 = 0.;
    let vy_max: f32 = 0.;

    // Init RNG, colors, forces and particles
    let mut rng: ThreadRng = rand::thread_rng();
    let colors: Vec<Color> = vec![
        Color::RED,
        Color::BLUE,
        Color::GREEN,
        Color::ORANGE,
        Color::YELLOW,
        Color::PURPLE,
        Color::BROWN,
        Color::PINK,
    ];
    let force_matrix = generate_force_matrix(colors.len(), &mut rng);
    let mut particles = generate_particles(NUM_PARTICLES, &mut rng, colors, x_min, x_max, y_min, y_max, vx_min, vx_max, vy_min, vy_max);

    // Init raylib
    let (mut rl, thread) = raylib::init()
        .size(params.width, params.height)
        .title("Particle Life")
        .build();

    // Drawing loop
    let mut current_time: f32 = 0.;
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);
        d.clear_background(Color::BLACK);
        // Update particles
        particles = update_particles(&particles, &params, &force_matrix);
        // // Mouse events
        // FIX!
        // if d.is_mouse_button_down(MouseButton::MOUSE_LEFT_BUTTON) {
        //     // Pick up particles
        //     let mut mouse_position = d.get_mouse_position();
        //     mouse_position.x /= params.width as f32;
        //     mouse_position.y /= params.height as f32;
        //     particles = particles
        //         .iter_mut()
        //         .map(|p| {
        //             p.velocity -= (p.position - mouse_position) * 0.5;
        //             *p
        //         })
        //         .collect();
        // }
        // Draw
        for p in &particles {
            d.draw_circle(p.position.x as i32, p.position.y as i32, p.mass, p.color);
        }
        d.draw_fps(4, 4);
        let text = format!("t = {:.2}", current_time);
        d.draw_text(&text, 4, 24, 18, Color::GRAY);
        current_time += params.time_step;
    }

}


/*
 TODO:
 - create submodules
 - gray outline on circles
 - Underlying grid
 - Space partitioning
 - Vectorize input conditions and boundaries
 - Change to seedable RNG to save a good state matrix
 - Add collisions
 - Color changes
 - Continuous color force
 - Background color potential
*/