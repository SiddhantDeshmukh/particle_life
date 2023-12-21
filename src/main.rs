use std::ffi::CStr;

use particle_life::*;

use rand::rngs::ThreadRng;
use raylib::{prelude::*, ffi::{GetMouseWheelMove, Rectangle}};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};


fn init(rng: &mut ThreadRng,
        window_width: i32, window_height: i32) -> Params {
    // Initialises the required Params for initial simulation based on
    // simulation mode
    let friction_half_life: f32 = 0.04;
    let time_step: f32 = 0.01;
    let max_radius: f32 = 0.1; // between 0 and 1
    let boundary_condition: BoundaryCondition = BoundaryCondition::Periodic;
    let grid_prgb_arr: [f32; 3] = generate_prgb_matrix(rng);
    let force_scale: f32 = 1.;

    return Params {
        window_width,
        window_height,
        time_step,
        friction_half_life,
        max_radius,
        boundary_condition,
        grid_prgb_arr,
        force_scale
    }
}

pub fn reset(rng: &mut ThreadRng, colors: &Vec<Color>, window_width: i32, window_height: i32) -> (Params, Vec<Particle>) {
    // Reinit, redraw
    let params = init(rng, window_width, window_height);
    let particles = generate_particles_cm(2000, rng, &colors, 0., 1., 0., 1.,
        0., 0., 0., 0.);
    (params, particles)
}

fn create_cstr(s: &str) -> Option<&CStr> {
    CStr::from_bytes_with_nul(s.as_bytes())
        .ok()
}

fn btn_rectangle(x: f32, y: f32) -> Rectangle {
    // Standard button rectangle
    Rectangle {x, y, width: 64., height: 32.}
}

fn slider(d: &mut RaylibDrawHandle, x: f32, y: f32, value: f32, min_value: f32, max_value: f32) -> f32 {
    // Wrapper around gui_slider()
    d.gui_slider(Rectangle {x, y, width: 96., height: 24.},
        create_cstr(&(min_value.to_string() + "\0")),
        create_cstr(&(max_value.to_string() + "\0")), value, min_value, max_value)
}

fn label(d: &mut RaylibDrawHandle, x: f32, y: f32, s: &str) {
    d.gui_label(Rectangle {x, y, width: 64., height: 24.}, create_cstr(s))
}

fn main() {
    // For rendering the window
    const WINDOW_WIDTH: i32 = 1600;
    const WINDOW_HEIGHT: i32 = 800;
    // Padding for aesthetics and control panel
    let window_left_padding = WINDOW_WIDTH / 6;  // holds control panel
    let window_right_padding = WINDOW_WIDTH / 48;  // aesthetic
    let window_top_padding = WINDOW_HEIGHT / 48;  // aesthetic
    let window_bot_padding = WINDOW_HEIGHT / 48;  // aesthetic
    let window_width_padding = window_left_padding + window_right_padding;
    let window_height_padding = window_top_padding + window_bot_padding;

    let sim_window_width: i32 = WINDOW_WIDTH - window_width_padding;
    let sim_window_height: i32 = WINDOW_HEIGHT - window_height_padding;

    // Init RNG
    let mut rng: ThreadRng = rand::thread_rng();
    // Simulation parameters
    let mut params = init(&mut rng, sim_window_width, sim_window_height);

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
    let colors: Vec<Color> = vec![
        Color::RED,
        Color::BLUE,
        Color::GREEN,
        Color::ORANGE,
        Color::YELLOW,
        Color::PURPLE,
        Color::BROWN,
        Color::PINK,
        Color::SKYBLUE,
        Color::LIME
    ];
    let force_matrix = generate_force_matrix(colors.len(), &mut rng);
    let mut particles = generate_particles_cm(NUM_PARTICLES, &mut rng, &colors, x_min, x_max, y_min, y_max, vx_min, vx_max, vy_min, vy_max);

    // Init raylib
    let (mut rl, thread) = raylib::init()
        .size(WINDOW_WIDTH, WINDOW_HEIGHT)
        .title("Particle Life")
        .build();

    let mut mouse_pickup_radius: f32 = 100.;

    // Simulation and drawing loop
    let mut current_time: f32 = 0.;
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);
        d.clear_background(Color::BLACK);
        // Update particles
        particles = update_particles_cm(&particles, &params, &force_matrix);
        // Drawing
        // Render pick-up circle around mouse
        let render_mouse_position = d.get_mouse_position();
        let mouse_position = Vector2 {x: render_mouse_position.x - window_left_padding as f32,
                                               y: render_mouse_position.y - window_top_padding as f32};
        // Pickup radius
        d.draw_circle_lines(render_mouse_position.x as i32, render_mouse_position.y as i32,
            mouse_pickup_radius, Color::GRAY);
        let world_mouse_position = win_to_world(mouse_position, params.window_width, params.window_height);
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
            // TODO:
            // Fix padding offsets
            particles = particles
                .par_iter_mut()
                .map(|p: &mut Particle| {
                    if mouse_position.distance_to(world_to_win(p.position, params.window_width, params.window_height)) <= mouse_pickup_radius {
                        p.velocity -= (p.position - world_mouse_position).normalized() * 0.5;
                    }
                    *p
                })
                .collect();
        }
        // Particles
        for p in &particles {
            d.draw_circle(p.win_pos_x(params.window_width) + window_left_padding,
                        p.win_pos_y(params.window_height) + window_top_padding,
                        3.,
                         p.color);
        }
        // Control panel
        // d.draw_rectangle(16, window_top_padding,
        //     window_left_padding - window_right_padding, WINDOW_HEIGHT - window_height_padding,
        //     Color::DARKGRAY);
        let text_x = 8.;
        // FPS
        d.draw_fps(text_x as i32, 18);
        // Timer
        let text = format!("t = {:.2}", current_time);
        d.draw_text(&text, text_x as i32, 40, 18, Color::WHITE);
        // Reset button
        let reset_button_clicked = d.gui_button(btn_rectangle(text_x as f32, 64.), create_cstr("Reset\0"));
        if reset_button_clicked {
            (params, particles) = reset(&mut rng, &colors, sim_window_width, sim_window_height);
            current_time = 0.;
        };
        // Params control
        let slider_x = text_x + 124.;
        label(&mut d, text_x as f32, 108., format!("Time step = {:.3}\0", params.time_step).as_str());
        params.time_step = slider(&mut d, slider_x, 108., params.time_step, 0.01, 0.1);
        label(&mut d, text_x as f32, 136., format!("Fric. 1/2-life = {:.2}\0", params.friction_half_life).as_str());
        params.friction_half_life = slider(&mut d, slider_x, 136., params.friction_half_life, 0.01, 0.99);
        label(&mut d, text_x as f32, 164., format!("Max Radius = {:.2}\0", params.max_radius).as_str());
        params.max_radius = slider(&mut d, slider_x, 164., params.max_radius, 0.01, 0.99);
        label(&mut d, text_x as f32, 192., format!("Force scale = {:.2}\0", params.force_scale).as_str());
        params.force_scale = slider(&mut d, slider_x, 192., params.force_scale, 0.1, 10.);
        // Boundary Condition
        // Force matrix

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
*/