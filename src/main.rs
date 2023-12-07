use particle_life::*;

use rand::rngs::ThreadRng;
use raylib::{prelude::*, ffi::GetMouseWheelMove};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};


fn init(rng: &mut ThreadRng,
        window_width: i32, window_height: i32,
        simulation_width: i32, simulation_height: i32) -> Params {
    // Initialises the required Params for initial simulation based on
    // simulation mode
    let friction_half_life: f32 = 0.5;
    let time_step: f32 = 0.05;
    let max_radius: f32 = 300.;
    let boundary_condition: BoundaryCondition = BoundaryCondition::Periodic;
    // let boundary_condition: BoundaryCondition = BoundaryCondition::Reflecting;
    let force_scale: f32 = 1.;
    let p1_rand_arr: [f32; 3] = generate_prgb_matrix(rng);
    let p2_rand_arr: [f32; 3] = generate_prgb_matrix(rng);

    return Params {
        window_width,
        window_height,
        simulation_width,
        simulation_height,
        friction_half_life,
        time_step,
        max_radius,
        boundary_condition,
        p1_rand_arr,
        p2_rand_arr,
        force_scale
    }
}

fn main() {
    // For rendering the window
    const WINDOW_WIDTH: i32 = 1600;
    const WINDOW_HEIGHT: i32 = 800;
    // Physics domain (for now the same size, need to decouple world from
    // drawing)
    let simulation_width: i32 = WINDOW_WIDTH;
    let simulation_height: i32 = WINDOW_HEIGHT;
    // Simulation parameters
    const NUM_PARTICLES: usize = 1000;  // TODO: make variable
    const MIN_MOUSE_PICKUP_RADIUS: f32 = 25.;
    const MAX_MOUSE_PICKUP_RADIUS: f32 = 500.;
    // const SEED: u64 = 420;
    // Init RNG and initial params
    let mut rng: ThreadRng = rand::thread_rng();
    let params = init(&mut rng, WINDOW_WIDTH, WINDOW_HEIGHT,
        simulation_width, simulation_height);

    // Bounds for initial particle spawning
    let x_min: f32 = params.x_min();
    let x_max: f32 = params.x_max();
    let y_min: f32 = params.y_min();
    let y_max: f32 = params.y_max();
    let vx_min: f32 = 0.;
    let vx_max: f32 = 0.;
    let vy_min: f32 = 0.;
    let vy_max: f32 = 0.;

    let mut particles = generate_particles(NUM_PARTICLES, &mut rng, x_min, x_max, y_min, y_max, vx_min, vx_max, vy_min, vy_max);

    // Init raylib
    let (mut rl, thread) = raylib::init()
        .size(params.window_width, params.window_height)
        .title("Particle Life")
        .build();

    let mut mouse_pickup_radius: f32 = 100.;

    // Simulatiuon and drawing loop
    let mut current_time: f32 = 0.;
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);
        d.clear_background(Color::BLACK);
        // Update particles
        particles = update_particles(&particles, &params);

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
            d.draw_circle(p.position.x as i32, p.position.y as i32, 4., p.color);
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
 - make RGB matrix a [f32; 9] 
 - Sim domain vs viewing domain
 - create submodules
 - Space partitioning
 - Vectorize input conditions and boundaries
 - Change to seedable RNG to save a good state matrix
 - Add collisions
 - Add spontaneous splitting
 - Continuous color force
    - Create 2 3-arrays (p1, p2) representing RGB
      Then compare
*/