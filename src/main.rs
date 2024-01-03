// use std::ffi::CStr;

use macroquad::ui::widgets::Checkbox;
use particle_life::*;

use macroquad::prelude::*;
use macroquad::prelude::Color;
use macroquad::prelude::MouseButton;

use ::rand::rngs::ThreadRng;
use raylib::{prelude::*, ffi::{GetMouseWheelMove, Rectangle}};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};


fn init(rng: &mut ThreadRng,
        window_width: f32, window_height: f32) -> Params {
    // Initialises the required Params for initial simulation based on
    // simulation mode
    let friction_half_life: f32 = 10.;
    let time_step: f32 = 5.;
    let max_radius: f32 = 20.; // between 0 and 100
    let boundary_condition: BoundaryCondition = BoundaryCondition::Periodic;
    let grid_prgb_arr: [f32; 3] = generate_prgb_matrix(rng);
    let force_scale: f32 = 3.;  // for user control

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

pub fn reset(rng: &mut ThreadRng, colors: &Vec<Color>, window_width: f32, window_height: f32) -> (Params, Vec<Particle>, Vec<Vec<f32>>) {
    // Reinit, redraw
    let params = init(rng, window_width, window_height);
    let particles = generate_particles_cm(2000, rng, &colors, 0., X_SIZE, 0., Y_SIZE,
        0., 0., 0., 0.);
    let force_matrix = generate_force_matrix(colors.len(), rng);
    (params, particles, force_matrix)
}

// fn create_cstr(s: &str) -> Option<&CStr> {
//     CStr::from_bytes_with_nul(s.as_bytes())
//         .ok()
// }

// fn btn_rectangle(x: f32, y: f32) -> Rectangle {
//     // Standard button rectangle
//     Rectangle {x, y, width: 64., height: 32.}
// }

// fn slider(d: &mut RaylibDrawHandle, x: f32, y: f32, value: f32, min_value: f32, max_value: f32) -> f32 {
//     // Wrapper around gui_slider()
//     d.gui_slider(Rectangle {x, y, width: 96., height: 24.},
//         create_cstr(&(min_value.to_string() + "\0")),
//         create_cstr(&(max_value.to_string() + "\0")), value, min_value, max_value)
// }

// fn label(d: &mut RaylibDrawHandle, x: f32, y: f32, s: &str) {
//     d.gui_label(Rectangle {x, y, width: 64., height: 24.}, create_cstr(s))
// }

// fn checkbox(d: &mut RaylibDrawHandle, x: f32, y: f32, checked: bool) -> bool {
//     d.gui_check_box(Rectangle {x, y, width: 24., height: 24.}, None, checked)
// }

// fn textbox(d: &mut RaylibDrawHandle, x: f32, y: f32, s: &str) -> i32 {
//     let mut text = Vec::from(s);
//     d.gui_text_input_box(Rectangle {x, y, width: 64., height: 24.}, None, None, None, &mut text)
// }

// fn scrollbar(d: &mut RaylibDrawHandle, x: f32, y: f32, width: f32, height: f32, value: i32, min_value: i32, max_value: i32) -> i32 {
//     d.gui_scroll_bar(Rectangle {x, y, width, height}, value, min_value, max_value)
// }

// fn get_first_char(s: &String) -> char {
//     s.chars().next().unwrap_or('\0')
// }

fn draw_fps(x: f32, y: f32, font_size: f32) {
    let fps = get_fps();
    let c = match fps {
        0..=10 => RED,
        11..=30 => ORANGE,
        _ => GREEN
    };
    draw_text(format!("FPS: {}", fps).as_str(), x, y, font_size, c)
}

#[macroquad::main("")]
async fn main() {
    // For rendering the window
    // const WINDOW_WIDTH: i32 = 1600;
    // const WINDOW_HEIGHT: i32 = 800;
    let mut window_width = screen_width();
    let mut window_height = screen_height();
    let mut seed = 420;  // for seedable RNG
    // Padding for aesthetics and control panel
    let window_left_padding = window_width / 6.;  // holds control panel
    let window_right_padding = window_width / 48.;  // aesthetic
    let window_top_padding = window_height / 48.;  // aesthetic
    let window_bot_padding = window_height / 48.;  // aesthetic
    let window_width_padding = window_left_padding + window_right_padding;
    let window_height_padding = window_top_padding + window_bot_padding;

    let sim_window_width: f32 = window_width - window_width_padding;
    let sim_window_height: f32 = window_height - window_height_padding;

    // Init RNG
    let mut rng: ThreadRng = ::rand::thread_rng();
    // let mut seed_rng = rand::SeedableRng::from_seed(seed);
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
    let colors: Vec<macroquad::color::Color> = vec![
        // Color::RED,
        // Color::ORANGE,
        // Color::YELLOW,
        // Color::LIME,
        // Color::PURPLE,
        // Color::BLUE,
        // Color::PINK,
        // Color::SKYBLUE,
        RED,
        ORANGE,
        YELLOW,
        LIME,
        PURPLE,
        BLUE,
        PINK,
        SKYBLUE
    ];
    // let color_names: Vec<String> = vec![
    //     String::from("Red"),
    //     String::from("Orange"),
    //     String::from("Yellow"),
    //     String::from("Lime"),
    //     String::from("Purple"),
    //     String::from("Blue"),
    //     String::from("Pink"),
    //     String::from("Skyblue"),
    // ];
    let mut force_matrix = generate_force_matrix(colors.len(), &mut rng);
    let mut particles = generate_particles_cm(NUM_PARTICLES, &mut rng, &colors, x_min, x_max, y_min, y_max, vx_min, vx_max, vy_min, vy_max);

    let mut mouse_pickup_radius: f32 = 100.;
    let mut is_paused: bool = false;
    let mut is_reversed: bool = false;
    let mut use_reflecting_bc: bool = false;  // by default use periodic

    let mut current_time: f32 = 0.;

    // Macroquad drawing loop
    loop {
        // Dynamic screen sizing
        window_width = screen_width();
        window_height = screen_height();
        // Drawing
        clear_background(BLACK);
        let text_x = 8.;
        // FPS
        draw_fps(text_x, 18., 32.);
        // Timer
        let text = format!("t = {:.2}", current_time);
        draw_text(text.as_str(), text_x, 32., 18., WHITE);
        // Reset button (resets with currently selected options)
        // let reset_button_clicked = d.gui_button(btn_rectangle(text_x + 132., 24.), create_cstr("Reset\0"));
        // if reset_button_clicked {
        //     (params, particles, force_matrix) = reset(&mut rng, &colors, sim_window_width, sim_window_height);
        //     current_time = 0.;
        // };
        // Pause button
        
        // label(&mut d, text_x, 64., "Paused: \0");
        // is_paused = checkbox(&mut d, text_x + 48., 64., is_paused);
        // // Reverse button
        // label(&mut d, text_x + 92., 64., "Reversed:\0");
        // is_reversed = checkbox(&mut d, text_x + 154., 64., is_reversed);
    //     // Params control
    //     let slider_x = text_x + 124.;
    //     label(&mut d, text_x, 108., format!("Time step = {:}\0", params.time_step).as_str());
    //     params.time_step = slider(&mut d, slider_x, 108., params.time_step, 1., 10.);
    //     label(&mut d, text_x, 136., format!("Fric. 1/2-life = {:}\0", params.friction_half_life).as_str());
    //     params.friction_half_life = slider(&mut d, slider_x, 136., params.friction_half_life, 1., 50.);
    //     label(&mut d, text_x, 164., format!("Max Radius = {:}\0", params.max_radius).as_str());
    //     params.max_radius = slider(&mut d, slider_x, 164., params.max_radius, 1., 100.);
    //     label(&mut d, text_x, 192., format!("Force scale = {:}\0", params.force_scale).as_str());
    //     params.force_scale = slider(&mut d, slider_x, 192., params.force_scale, 1., 10.);
    //     // Boundary Condition
    //     label(&mut d, text_x, 220., "Reflecting BC? (otherwise periodic)\0");
    //     use_reflecting_bc = checkbox(&mut d, slider_x + 80., 220., use_reflecting_bc);
    //     if use_reflecting_bc {
    //         params.boundary_condition = BoundaryCondition::Reflecting;
    //     } else {
    //         params.boundary_condition = BoundaryCondition::Periodic;
    //     }

    //     // Force matrix
    //     // TODO: Scrollbox, one row for each color-color interaction
    //     let base_y: f32 = 240.;
    //     // Color reference
    //     let cref_text = "R=Red, O=Orange, Y=Yellow, L=Lime, P=Purple,\nB=Blue, P=Pink, S=Skyblue\0";
    //     label(&mut d, text_x, base_y, &cref_text);

    //     let mut y_increment: f32 = 32.;
    //     for ci in 0..color_names.len() {
    //         for cj in 0..color_names.len() {
    //             let text = format!("{:?} - {:?}: {:.2}\0", get_first_char(&color_names[ci]), get_first_char(&color_names[cj]), force_matrix[ci][cj]);
    //             // TODO: Set colors and add sliders
    //             label(&mut d, text_x, base_y + y_increment, text.as_str());
    //             force_matrix[ci][cj] = slider(&mut d, text_x + 128., base_y + y_increment, force_matrix[ci][cj], -1., 1.);
    //             y_increment += 12.
    //         }
    //     }

    //     // Render pick-up circle around mouse
        let render_mouse_position = Vec2::from(mouse_position());
        let mouse_position_v = Vector2 {x: render_mouse_position.x - window_left_padding as f32,
                                               y: render_mouse_position.y - window_top_padding as f32};
        // Pickup radius
        draw_circle_lines(render_mouse_position.x, render_mouse_position.y, mouse_pickup_radius, 1., GRAY);
        let world_mouse_position = win_to_world(mouse_position_v, params.window_width, params.window_height);
    //     let mouse_wheel = unsafe { GetMouseWheelMove() };
    //     // Mouse events
    //     // Mouse wheel to change pickup size
        let mouse_wheel_y = mouse_wheel().1;
        if mouse_wheel_y != 0. {
            mouse_pickup_radius = clamp(mouse_pickup_radius + (mouse_wheel_y * 25.),
            MIN_MOUSE_PICKUP_RADIUS, MAX_MOUSE_PICKUP_RADIUS)
        }
        // Particles
        // No physics updates when paused
        if !is_paused {
            // Click to grab
            // TODO:
            // Handle this logic in the update loop to prevent double looping
            if is_mouse_button_down(MouseButton::Left) {
                // Pick up particles near mouse
                particles = particles
                    .par_iter_mut()
                    .map(|p: &mut Particle| {
                        
                        if mouse_position_v.distance_to(world_to_win(p.position, params.window_width, params.window_height)) <= mouse_pickup_radius {
                            p.velocity -= (p.position - world_mouse_position).normalized() * 0.5;
                        }
                        *p
                    })
                    .collect();
            }
            particles = update_particles_cm(&particles, &params, &force_matrix, false);
            if is_reversed {
                current_time -= params.time_step;
            } else {
                current_time += params.time_step;
            }
        }
        // Draw updated particles
        for p in &particles {
            draw_circle((p.win_pos_x(window_width) + window_left_padding) as f32,
            (p.win_pos_y(window_height) + window_top_padding) as f32,
            3., p.color)
        }
        // for p in &particles {
        //     d.draw_circle(p.win_pos_x(params.window_width) + window_left_padding,
        //                 p.win_pos_y(params.window_height) + window_top_padding,
        //                 3.,
        //                  p.color);
        // }

        // Next frame
        next_frame().await
    }
    // // Init raylib
    // let (mut rl, thread) = raylib::init()
    //     .size(WINDOW_WIDTH, WINDOW_HEIGHT)
    //     .title("Particle Life")
    //     .build();
    // while !rl.window_should_close() {
    //     // Drawing
    //     let mut d = rl.begin_drawing(&thread);
    //     d.clear_background(Color::BLACK);
    //     let text_x = 8.;
    //     // Draw UI
    //     // Control panel
    //     // d.draw_rectangle(16, window_top_padding,
    //     //     window_left_padding - window_right_padding, WINDOW_HEIGHT - window_height_padding,
    //     //     Color::DARKGRAY);
    //     // FPS
    //     d.draw_fps(text_x as i32, 18);
    //     // Timer
    //     let text = format!("t = {:}", current_time);
    //     d.draw_text(&text, text_x as i32, 40, 18, Color::WHITE);
    //     // Pause button
    //     label(&mut d, text_x, 64., "Paused: \0");
    //     is_paused = checkbox(&mut d, text_x + 48., 64., is_paused);
    //     // Reverse button
    //     label(&mut d, text_x + 92., 64., "Reversed:\0");
    //     is_reversed = checkbox(&mut d, text_x + 154., 64., is_reversed);
    //     // Params control
    //     let slider_x = text_x + 124.;
    //     label(&mut d, text_x, 108., format!("Time step = {:}\0", params.time_step).as_str());
    //     params.time_step = slider(&mut d, slider_x, 108., params.time_step, 1., 10.);
    //     label(&mut d, text_x, 136., format!("Fric. 1/2-life = {:}\0", params.friction_half_life).as_str());
    //     params.friction_half_life = slider(&mut d, slider_x, 136., params.friction_half_life, 1., 50.);
    //     label(&mut d, text_x, 164., format!("Max Radius = {:}\0", params.max_radius).as_str());
    //     params.max_radius = slider(&mut d, slider_x, 164., params.max_radius, 1., 100.);
    //     label(&mut d, text_x, 192., format!("Force scale = {:}\0", params.force_scale).as_str());
    //     params.force_scale = slider(&mut d, slider_x, 192., params.force_scale, 1., 10.);
    //     // Boundary Condition
    //     label(&mut d, text_x, 220., "Reflecting BC? (otherwise periodic)\0");
    //     use_reflecting_bc = checkbox(&mut d, slider_x + 80., 220., use_reflecting_bc);
    //     if use_reflecting_bc {
    //         params.boundary_condition = BoundaryCondition::Reflecting;
    //     } else {
    //         params.boundary_condition = BoundaryCondition::Periodic;
    //     }

    //     // Force matrix
    //     // TODO: Scrollbox, one row for each color-color interaction
    //     let base_y: f32 = 240.;
    //     // Color reference
    //     let cref_text = "R=Red, O=Orange, Y=Yellow, L=Lime, P=Purple,\nB=Blue, P=Pink, S=Skyblue\0";
    //     label(&mut d, text_x, base_y, &cref_text);

    //     let mut y_increment: f32 = 32.;
    //     for ci in 0..color_names.len() {
    //         for cj in 0..color_names.len() {
    //             let text = format!("{:?} - {:?}: {:.2}\0", get_first_char(&color_names[ci]), get_first_char(&color_names[cj]), force_matrix[ci][cj]);
    //             // TODO: Set colors and add sliders
    //             label(&mut d, text_x, base_y + y_increment, text.as_str());
    //             force_matrix[ci][cj] = slider(&mut d, text_x + 128., base_y + y_increment, force_matrix[ci][cj], -1., 1.);
    //             y_increment += 12.
    //         }
    //     }

    //     // Render pick-up circle around mouse
    //     let render_mouse_position = d.get_mouse_position();
    //     let mouse_position = Vector2 {x: render_mouse_position.x - window_left_padding as f32,
    //                                            y: render_mouse_position.y - window_top_padding as f32};
    //     // Pickup radius
    //     d.draw_circle_lines(render_mouse_position.x as i32, render_mouse_position.y as i32,
    //         mouse_pickup_radius, Color::GRAY);
    //     let world_mouse_position = win_to_world(mouse_position, params.window_width, params.window_height);
    //     let mouse_wheel = unsafe { GetMouseWheelMove() };
    //     // Mouse events
    //     // Mouse wheel to change pickup size
    //     if mouse_wheel != 0. {
    //         mouse_pickup_radius = clamp(mouse_pickup_radius + (mouse_wheel * 25.),
    //         MIN_MOUSE_PICKUP_RADIUS, MAX_MOUSE_PICKUP_RADIUS)
    //     }
    //     // Particles
    //     // No physics updates when paused
    //     if !is_paused {
    //         // Click to grab
    //         // TODO:
    //         // Handle this logic in the update loop to prevent double looping
    //         if d.is_mouse_button_down(MouseButton::MOUSE_LEFT_BUTTON) {
    //             // Pick up particles near mouse
    //             particles = particles
    //                 .par_iter_mut()
    //                 .map(|p: &mut Particle| {
    //                     if mouse_position.distance_to(world_to_win(p.position, params.window_width, params.window_height)) <= mouse_pickup_radius {
    //                         p.velocity -= (p.position - world_mouse_position).normalized() * 0.5;
    //                     }
    //                     *p
    //                 })
    //                 .collect();
    //         }
    //         particles = update_particles_cm(&particles, &params, &force_matrix, is_reversed);
    //         if is_reversed {
    //             current_time -= params.time_step;
    //         } else {
    //             current_time += params.time_step;
    //         }
    //     }

    //     for p in &particles {
    //         d.draw_circle(p.win_pos_x(params.window_width) + window_left_padding,
    //                     p.win_pos_y(params.window_height) + window_top_padding,
    //                     3.,
    //                      p.color);
    //     }
    //     // Reset button (resets with currently selected options)
    //     let reset_button_clicked = d.gui_button(btn_rectangle(text_x + 132., 24.), create_cstr("Reset\0"));
    //     if reset_button_clicked {
    //         (params, particles, force_matrix) = reset(&mut rng, &colors, sim_window_width, sim_window_height);
    //         current_time = 0.;
    //     };
    // }

}


/*
 TODO:
 - render with macroquad instead
 - dynamic resize UI (draw UI func)
 - add/remove particles with brush
 - scrollbox for force matrix
 - create submodules
 - Space partitioning
 - Vectorize input conditions and boundaries
 - Change to seedable RNG to save a good state matrix, then reset with seed
 - Add collisions
 - Add spontaneous splitting
*/