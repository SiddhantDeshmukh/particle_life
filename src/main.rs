use particle_life::*;

use rand::rngs::ThreadRng;
use raylib::{prelude::*, ffi::GetMouseWheelMove};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

fn main() {
    // For rendering the window
    const WINDOW_WIDTH: i32 = 1600;
    const WINDOW_HEIGHT: i32 = 800;

    // Physics domain (for now the same size, need to decouple world from
    // drawing)
    let simulation_width: i32 = WINDOW_WIDTH;
    let simulation_height: i32 = WINDOW_HEIGHT;
    // Simulation parameters
    let params = Params {
        width: simulation_width,
        height: simulation_height,
        friction_half_life: 0.5,
        time_step: 0.1,
        // time_step: 0.5,  // for debugging
        max_radius: 200.,
        boundary_condition: BoundaryCondition::Periodic
        // boundary_condition: BoundaryCondition::Reflecting
    };

    // let window_padding = 100;

    const NUM_PARTICLES: usize = 2000;  // TODO: make variable
    const MIN_MOUSE_PICKUP_RADIUS: f32 = 25.;
    const MAX_MOUSE_PICKUP_RADIUS: f32 = 500.;
    // const SEED: u64 = 420;

    // Bounds for initial particle spawning
    let x_min: f32 = params.x_min();
    let x_max: f32 = params.x_max();
    let y_min: f32 = params.y_min();
    let y_max: f32 = params.y_max();
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
    let mut particles = generate_particles_cm(NUM_PARTICLES, &mut rng, colors, x_min, x_max, y_min, y_max, vx_min, vx_max, vy_min, vy_max);

    // Init raylib
    let (mut rl, thread) = raylib::init()
        .size(WINDOW_WIDTH, WINDOW_HEIGHT)
        .title("Particle Life")
        .build();

    let mut mouse_pickup_radius: f32 = 100.;

    // Simulatiuon and drawing loop
    let mut current_time: f32 = 0.;
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);
        d.clear_background(Color::BLACK);
        // Update particles
        particles = update_particles_cm(&particles, &params, &force_matrix);

        // Render pick-up circle around mouse
        let mouse_position = d.get_mouse_position();
        let mouse_wheel = unsafe { GetMouseWheelMove() };
        // Mouse events
        // Mouse wheel to change pickup size
        if mouse_wheel != 0. {
            mouse_pickup_radius = clamp(mouse_pickup_radius + (mouse_wheel * 25.),
            MIN_MOUSE_PICKUP_RADIUS, MAX_MOUSE_PICKUP_RADIUS)
        }
        // Click to grab
        if d.is_mouse_button_down(MouseButton::MOUSE_LEFT_BUTTON) {
            // Pick up particles near mouse
            particles = particles
                .par_iter_mut()
                .map(|p| {
                    if mouse_position.distance_to(p.position) <= mouse_pickup_radius {
                        p.velocity -= (p.position - mouse_position) * 0.5;
                    }
                    *p
                })
                .collect();
        }
        // Draw
        // Particles
        for p in &particles {
            d.draw_circle(p.position.x as i32, p.position.y as i32, 3., p.color);
        }
        // FPS
        d.draw_fps(4, 4);
        // Pickup radius
        d.draw_circle_lines(mouse_position.x as i32, mouse_position.y as i32,
            mouse_pickup_radius, Color::GRAY);
        // Timer
        let text = format!("t = {:.2}", current_time);
        d.draw_text(&text, 4, 24, 18, Color::GRAY);
        current_time += params.time_step;
    }

}


/*
 TODO:
 - User controlled params
    - friction, dt
    - change forces
 - add/remove particles with brush
 - Sim domain vs viewing domain
 - create submodules
 - Space partitioning
 - Vectorize input conditions and boundaries
 - Change to seedable RNG to save a good state matrix
 - Add collisions
 - Add spontaneous splitting
 - Continuous color force
*/