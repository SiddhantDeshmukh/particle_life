use egui_macroquad::egui::Widget;
use particle_life::*;

use macroquad::prelude::*;
use macroquad::prelude::Color;
use macroquad::prelude::MouseButton;

use egui_macroquad::egui;

use ::rand::rngs::ThreadRng;
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

fn draw_fps(x: f32, y: f32, font_size: f32) {
    let fps = get_fps();
    let c = match fps {
        0..=10 => RED,
        11..=30 => ORANGE,
        _ => GREEN
    };
    draw_text(format!("FPS: {}", fps).as_str(), x, y, font_size, c)
}

macro_rules! slider {
    ($value:expr, $min:expr, $max:expr, $ui:expr) => {
        egui::Slider::new($value, std::ops::RangeInclusive::new($min, $max)).ui($ui);
    };
}

#[macroquad::main("ParticleLife")]
async fn main() {
    // For rendering the window
    let mut window_width = screen_width();
    let mut window_height = screen_height();
    // let mut seed = 420;  // for seedable RNG

    let sim_window_width: f32 = window_width;
    let sim_window_height: f32 = window_height;
    let window_color = egui::Color32::from_rgba_unmultiplied(8, 0, 72, 200);

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
        RED,
        ORANGE,
        YELLOW,
        LIME,
        PURPLE,
        BLUE,
        PINK,
        SKYBLUE
    ];
    let mut force_matrix = generate_force_matrix(colors.len(), &mut rng);
    let mut particles = generate_particles_cm(NUM_PARTICLES, &mut rng, &colors, x_min, x_max, y_min, y_max, vx_min, vx_max, vy_min, vy_max);

    let mut mouse_pickup_radius: f32 = 100.;
    let mut is_paused: bool = false;
    let mut is_reversed: bool = false;

    let mut current_time: f32 = 0.;

    // Macroquad drawing loop with egui
    loop {
        // Dynamic screen sizing
        window_width = screen_width();
        window_height = screen_height();
        // Determine font size based on scren size
        let font_size_medium: f32 = 0.05 * window_height;
        // Drawing
        clear_background(BLACK);
        let text_x = 8.;

        // Render pick-up circle around mouse
        let render_mouse_position = Vec2::from(mouse_position());
        let mouse_position_v = Vec2::new(render_mouse_position.x, render_mouse_position.y);
        // Pickup radius
        draw_circle_lines(render_mouse_position.x, render_mouse_position.y, mouse_pickup_radius, 1., GRAY);
        let world_mouse_position = win_to_world(mouse_position_v, params.window_width, params.window_height);
        // Mouse events
        // Mouse wheel to change pickup size
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
                        if mouse_position_v.distance(world_to_win(p.position, params.window_width, params.window_height)) <= mouse_pickup_radius {
                            p.velocity -= (p.position - world_mouse_position).normalize() * 0.5;
                        }
                        *p
                    })
                    .collect();
            }
            particles = update_particles_cm(&particles, &params, &force_matrix, is_reversed);
            if is_reversed {
                current_time -= params.time_step;
            } else {
                current_time += params.time_step;
            }
        }
        // Draw updated particles
        for p in &particles {
            draw_circle(p.win_pos_x(window_width), p.win_pos_y(window_height),
            2.5, p.color)
        }
        // UI Window
        egui_macroquad::ui(|egui_ctx| {
            egui::Window::new("Controls")
                .frame(egui::Frame{fill:window_color, ..Default::default()})
                .fixed_pos(egui::pos2(text_x, 0.05*window_height))
                // .fixed_size(egui::vec2(window_width * 0.1, window_height * 0.4))
                .show(egui_ctx, |ui| {
                    egui_ctx.set_pixels_per_point(0.0015 * window_height);
                    ui.horizontal(|ui| {
                        egui::FontId::default().size = 0.02 * window_height;
                        // Reset button
                        if ui.button("Reset").clicked() {
                            (params, particles, force_matrix) = reset(&mut rng, &colors, sim_window_width, sim_window_height);
                        }
                        // Pause & reverse checkboxes
                        ui.checkbox(&mut is_paused, "Pause");
                        ui.checkbox(&mut is_reversed, "Reverse");
                    });
                    // Params control
                    ui.horizontal(|ui| {
                        ui.label("Time step");
                        slider!(&mut params.time_step, 1., 10., ui);
                    });
                    ui.horizontal(|ui| {
                        ui.label("Friction 1/2-life");
                        slider!(&mut params.friction_half_life, 1., 50., ui);
                    });
                    ui.horizontal(|ui| {
                        ui.label("Max Interaction Radius");
                        slider!(&mut params.max_radius, 1., 100., ui);
                    });
                    ui.horizontal(|ui| {
                        ui.label("Force scale");
                        slider!(&mut params.force_scale, 1., 10., ui);
                    });
                    // Boundary conditions
                    egui::ComboBox::from_label("Boundary Condition")
                        .selected_text(format!("{:?}", params.boundary_condition))
                        .show_ui(ui, |ui| {
                            // TODO:
                            // - add "pushing" condition that draws to centre
                            ui.selectable_value(&mut params.boundary_condition, BoundaryCondition::Periodic, "Periodic");
                            ui.selectable_value(&mut params.boundary_condition, BoundaryCondition::Reflecting, "Reflecting");
                    });
                    // Force matrix
                    ui.collapsing("Force Matrix", |ui| {
                        egui::ScrollArea::new([false, true]).show(ui, |ui| {
                            for ci in 0..colors.len() {
                                for cj in 0..colors.len() {
                                    ui.horizontal(|ui| {
                                        ui.colored_label(egui_color(colors[ci]), "⬜");
                                        ui.label("->");
                                        ui.colored_label(egui_color(colors[cj]), "⬜");
                                        slider!(&mut force_matrix[ci][cj], -1., 1., ui);
                                    });
                                }
                            }
                        })
                    });
                });
        });
        // Draw control panel
        egui_macroquad::draw();
        // FPS and Timer
        draw_fps(text_x, 0.03 * window_height, font_size_medium);
        let text = format!("t = {:.2}", current_time);
        draw_text(text.as_str(), text_x + 0.15*window_width, 0.03 * window_height,
            font_size_medium, WHITE);
        // Next frame
        next_frame().await
    }
}


/*
 TODO:
 - dynamic resize UI (draw UI func)
 - add/remove particles with brush
 - Save/Load config
 - create submodules
 - Space partitioning
 - Vectorize input conditions and boundaries
 - Change to seedable RNG to save a good state matrix, then reset with seed
 - Add collisions
 - Add spontaneous splitting
*/