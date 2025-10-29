// ============================================================================
// RVM (Rotating Vector Model) Visualization
// Visualizes pulsar magnetic field lines, polarization vectors, and
// polarization angle/intensity plots
// ============================================================================

// ============================================================================
// IMPORTS
// ============================================================================
import * as THREE from '../node_modules/three/build/three.module.js';
import { GUI } from '../node_modules/dat.gui/build/dat.gui.module.js';
import { OrbitControls } from '../node_modules/three/examples/jsm/controls/OrbitControls.js';
import { Line2 } from '../node_modules/three/examples/jsm/lines/Line2.js';
import { LineMaterial } from '../node_modules/three/examples/jsm/lines/LineMaterial.js';
import { LineGeometry } from '../node_modules/three/examples/jsm/lines/LineGeometry.js';

// ============================================================================
// CONFIGURATION
// ============================================================================

// Physical and visualization parameters
const Config = function () {
  // Physical parameters
  this.dipole_angle = 30.0;        // Angle between magnetic and rotation axes (degrees)
  this.obs_angle = 90.0;           // Observer's viewing angle from rotation axis (degrees)
  this.pl_radius = 10.0;           // Radius of polarization limiting sphere
  this.emission_h = 0.01;          // Height of radio emission region
  this.translation_d = 0.41;       // Offset distance of emission beam
  this.translation_delta = 0.11;   // Angular aperture of emission beam (degrees)
  
  // Animation parameters
  this.rotation_speed = 0.01;      // Speed of pulsar rotation
  
  // Display modes
  this.x_mode = false;             // Show X-mode (extraordinary) polarization
  this.o_mode = true;              // Show O-mode (ordinary) polarization
  this.rotating = false;           // Enable/disable rotation animation
  this.lock_observer = true;       // Lock observer viewpoint
};

const conf = new Config();

// Quaternions for dipole rotation transformations
const q_dipole = new THREE.Quaternion();
const q_dipole_inv = new THREE.Quaternion();

// Update quaternions based on current dipole angle
function updateQuaternion() {
  q_dipole.setFromAxisAngle(new THREE.Vector3(1, 0, 0), Math.PI * conf.dipole_angle / 180.0);
  q_dipole_inv.copy(q_dipole).invert();
}

updateQuaternion();

// ============================================================================
// THREE.JS SCENE SETUP
// ============================================================================

// Renderer setup
const width = window.innerWidth;
const height = window.innerHeight;
const canvas = document.getElementById("vis");
const renderer = new THREE.WebGLRenderer({
  canvas: canvas,
  alpha: true,
  antialias: true,
  preserveDrawingBuffer: true
});
renderer.setSize(width, height);

// Scene setup
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x000000);

// Camera setup
const camera = new THREE.PerspectiveCamera(30, width / height, 1, 1000);
const camera_distance = 160;
camera.position.z = camera_distance * Math.cos(80 * Math.PI / 180);
camera.position.x = camera_distance * Math.sin(80 * Math.PI / 180);
camera.up.set(0, 0, 1);
camera.lookAt(new THREE.Vector3(0, 0, 0));

// Camera controls
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.05;
controls.enableKeys = false;
controls.minPolarAngle = controls.maxPolarAngle = conf.obs_angle * Math.PI / 180;

// Lighting
const directionalLight = new THREE.DirectionalLight(0xffffff);
directionalLight.position.set(107, 107, 107);
scene.add(directionalLight);

const ambientLight = new THREE.AmbientLight(0x404040);
scene.add(ambientLight);

// ============================================================================
// 3D OBJECTS - Neutron Star and Field Visualization
// ============================================================================

// Neutron star sphere
const star_geometry = new THREE.SphereGeometry(1.0, 64, 64);
const star_material = new THREE.MeshPhongMaterial({ color: 0xaaaaaa });
const star_mesh = new THREE.Mesh(star_geometry, star_material);
scene.add(star_mesh);

// Polarization limiting sphere (translucent)
const pl_geometry = new THREE.SphereGeometry(1.0, 64, 64);
const pl_material = new THREE.MeshPhongMaterial({ 
  color: 0xffffff,
  transparent: true,
  depthWrite: false,
  blending: THREE.NormalBlending,
  opacity: 0.3
});
const pl_mesh = new THREE.Mesh(pl_geometry, pl_material);
pl_mesh.scale.setScalar(conf.pl_radius);
scene.add(pl_mesh);

// Rotation axis (vertical line)
const axis_geometry = new THREE.BufferGeometry();
const axis_vertices = new Float32Array([0, 0, -300, 0, 0, 300]);
axis_geometry.setAttribute('position', new THREE.BufferAttribute(axis_vertices, 3));
const spin_axis = new THREE.Line(axis_geometry, new THREE.LineBasicMaterial({
  color: 0x8080aa,
  linewidth: 2.5,
}));
scene.add(spin_axis);

// Groups for dynamic field lines and polarization vectors
const field_lines = new THREE.Group();
const polarization_vectors = new THREE.Group();
const intersect_points = []; // Stores intersection points for polarization calculations

scene.add(field_lines);
scene.add(polarization_vectors);

// ============================================================================
// EMISSION CONE VISUALIZATION
// ============================================================================

// Load texture for emission cone
const textureLoader = new THREE.TextureLoader();
const plasma_texture = textureLoader.load('plasma.jpg');

// Create emission cone visualization
function createEmissionCone() {
  const cone_material = new THREE.MeshPhongMaterial({
    map: plasma_texture,
    transparent: true,
    opacity: 0.5,
    depthWrite: false,
    side: THREE.DoubleSide,
    alphaTest: 0.1
  });

  // Compute cone position based on translation parameters
  const deltaRad = conf.translation_delta * Math.PI / 180;
  const px = conf.translation_d * Math.sin(deltaRad) * Math.cos(0);
  const py = conf.translation_d * Math.sin(deltaRad) * Math.sin(0);
  const pz = conf.translation_d * Math.cos(deltaRad);
  
  // Upper cone
  const cone_geometry_upper = new THREE.ConeGeometry(5, 100, 64, 32, true);
  const cone_upper = new THREE.Mesh(cone_geometry_upper, cone_material);
  cone_upper.name = "line";
  cone_upper.rotateX(-(conf.dipole_angle / 180 * Math.PI + Math.PI / 2));
  cone_upper.position.set(px, py, pz);
  cone_geometry_upper.translate(0, -50, 0);
  field_lines.add(cone_upper);

  // Lower cone
  const cone_geometry_lower = new THREE.ConeGeometry(5, 100, 50);
  const cone_lower = new THREE.Mesh(cone_geometry_lower, cone_material);
  cone_lower.name = "line";
  cone_lower.rotateX(-(conf.dipole_angle / 180 * Math.PI - Math.PI / 2));
  cone_lower.position.set(px, py, pz);
  cone_geometry_lower.translate(0, -50, 0);
  field_lines.add(cone_lower);
}

// ============================================================================
// OBSERVER TRAJECTORY CIRCLE
// ============================================================================

// Create circle showing observer's trajectory at polarization limiting radius
const obs_circle_positions = [];
const obs_circle_segments = 80;
for (let i = 0; i <= obs_circle_segments; i++) {
  const phi = i * 2.0 * Math.PI / obs_circle_segments;
  obs_circle_positions.push(Math.cos(phi), Math.sin(phi), 0.0);
}

const obs_circle_geometry = new LineGeometry();
obs_circle_geometry.setPositions(obs_circle_positions);

const obs_circle_material = new LineMaterial({
  color: 0x846204,
  worldUnits: true,
  linewidth: 0.10,
  vertexColors: false,
  dashed: false,
  transparent: true,
  opacity: 0.7,
  blending: THREE.CustomBlending,
});

const obs_circle = new Line2(obs_circle_geometry, obs_circle_material);
obs_circle.computeLineDistances();
scene.add(obs_circle);
obs_circle.scale.setScalar((conf.pl_radius + conf.emission_h) * Math.sin(conf.obs_angle * Math.PI / 180));
obs_circle.position.setZ((conf.pl_radius + conf.emission_h) * Math.cos(conf.obs_angle * Math.PI / 180));

// ============================================================================
// MAGNETIC FIELD LINES
// ============================================================================

// Clear all objects with specific name from a group
function clearGroup(group) {
  for (let i = group.children.length - 1; i >= 0; i--) {
    const obj = group.children[i];
    if (obj.name === "line" || obj.name === "arrow") {
      if (obj.material) obj.material.dispose();
      if (obj.geometry) obj.geometry.dispose();
      group.remove(obj);
    }
  }
}

// Create dipole magnetic field lines
function createFieldLines(theta_obs, radius_pl, num_lines, num_segments, offset_d, offset_delta) {
  intersect_points.length = 0;
  
  for (let i = 0; i < num_lines; i++) {
    const phi = i * 2.0 * Math.PI / num_lines;
    
    // Calculate offset vector
    const offset = new THREE.Vector3(
      offset_d * Math.sin(offset_delta) * Math.cos(phi),
      offset_d * Math.sin(offset_delta) * Math.sin(phi),
      offset_d * Math.cos(offset_delta)
    );

    // Intersection point at polarization limiting radius (lab frame)
    const p_intersect_lab = new THREE.Vector3(
      radius_pl * Math.sin(theta_obs) * Math.cos(phi),
      radius_pl * Math.sin(theta_obs) * Math.sin(phi),
      radius_pl * Math.cos(theta_obs)
    );
    
    intersect_points.push(p_intersect_lab.clone());
    
    // Transform to dipole frame
    const p_intersect = p_intersect_lab.clone().sub(offset);
    p_intersect.applyQuaternion(q_dipole);
    
    // Calculate field line parameters
    const r_local = p_intersect.length();
    const cosArg = Math.max(-1, Math.min(1, p_intersect.z / r_local));
    const theta_intersect = Math.acos(cosArg);
    
    // Calculate L-shell parameter
    const L = r_local / (Math.sin(theta_intersect) ** 2);
    const phi_line = Math.atan2(p_intersect.y, p_intersect.x);
    
    // Generate field line points
    const line_positions = [];
    for (let j = 0; j <= num_segments; j++) {
      const theta = j * Math.PI / num_segments;
      const r = L * Math.sin(theta) ** 2;
      const p = new THREE.Vector3(
        r * Math.sin(theta) * Math.cos(phi_line),
        r * Math.sin(theta) * Math.sin(phi_line),
        r * Math.cos(theta)
      );
      // Transform back to lab frame
      p.applyQuaternion(q_dipole_inv);
      p.add(offset);
      line_positions.push(p.x, p.y, p.z);
    }
    
    // Create field line geometry
    const line_geometry = new LineGeometry();
    line_geometry.setPositions(line_positions);
    
    // Field line material
    const line_material = new LineMaterial({
      color: 0x024A0D,
      vertexColors: false,
      worldUnits: true,
      linewidth: 0.10,
      transparent: true,
      opacity: 0.7,
      alphaToCoverage: true,
      blending: THREE.AdditiveBlending,
      blendEquation: THREE.AddEquation,
      blendSrc: THREE.SrcAlphaFactor,
      blendDst: THREE.OneMinusSrcAlphaFactor,
    });

    const field_line = new Line2(line_geometry, line_material);
    field_line.computeLineDistances();
    field_line.name = "line";
    field_lines.add(field_line);
  }
  
  createEmissionCone();
}

// ============================================================================
// POLARIZATION VECTORS
// ============================================================================

// Create 3D polarization arrows at intersection points
function createPolarizationVectors() {
  const dipole_moment = new THREE.Vector3(
    0,
    Math.sin(conf.dipole_angle * Math.PI / 180),
    Math.cos(conf.dipole_angle * Math.PI / 180)
  );
  const radius_cubed = Math.pow(conf.pl_radius * 10, 3);
  const radius_fifth = Math.pow(conf.pl_radius * 10, 5);

  for (const p of intersect_points) {
    // Calculate magnetic field direction
    const b = p.clone().multiplyScalar(3.0 * dipole_moment.dot(p) / radius_fifth);
    const b2 = dipole_moment.clone().multiplyScalar(1.0 / radius_cubed);
    b.sub(b2);
    
    // Calculate polarization direction
    let dir = new THREE.Vector3();
    dir.crossVectors(p, b);
    dir.normalize();
    
    if (conf.o_mode) {
      dir.cross(p);
      dir.normalize();
    }

    // Create arrow geometry (shaft + cone head)
    const arrow_length = 3.0;
    const head_length = Math.min(0.4 * arrow_length, 0.5);
    const shaft_length = arrow_length - head_length;
    const shaft_radius = 0.05;
    const head_radius = 0.1;

    const arrow_material = new THREE.MeshPhongMaterial({ color: 0xBF3A0A });

    // Shaft (cylinder)
    const shaft_geometry = new THREE.CylinderGeometry(shaft_radius, shaft_radius, shaft_length, 8);
    shaft_geometry.translate(0, shaft_length / 2, 0);
    const shaft_mesh = new THREE.Mesh(shaft_geometry, arrow_material);

    // Head (cone)
    const head_geometry = new THREE.ConeGeometry(head_radius, head_length, 10);
    head_geometry.translate(0, shaft_length + head_length / 2, 0);
    const head_mesh = new THREE.Mesh(head_geometry, arrow_material);

    // Combine into arrow group
    const arrow = new THREE.Group();
    arrow.add(shaft_mesh);
    arrow.add(head_mesh);

    // Orient arrow along polarization direction
    const up = new THREE.Vector3(0, 1, 0);
    const quaternion = new THREE.Quaternion().setFromUnitVectors(up, dir.clone().normalize());
    arrow.setRotationFromQuaternion(quaternion);

    arrow.position.copy(p);
    arrow.name = "arrow";
    polarization_vectors.add(arrow);
  }
}

// ============================================================================
// UPDATE FUNCTIONS
// ============================================================================

// Update polarization limiting sphere and field lines
function updatePolarizationSphere() {
  pl_mesh.scale.setScalar(conf.pl_radius);
  obs_circle.scale.setScalar((conf.pl_radius + conf.emission_h) * Math.sin(conf.obs_angle * Math.PI / 180));
  obs_circle.position.setZ((conf.pl_radius + conf.emission_h) * Math.cos(conf.obs_angle * Math.PI / 180));
  controls.minPolarAngle = controls.maxPolarAngle = conf.obs_angle * Math.PI / 180;
  updateFieldLines();
}

// Regenerate field lines and polarization vectors
function updateFieldLines() {
  updateQuaternion();
  clearGroup(field_lines);
  clearGroup(polarization_vectors);
  createFieldLines(
    conf.obs_angle * Math.PI / 180,
    conf.pl_radius + conf.emission_h,
    40,
    100,
    conf.translation_d,
    conf.translation_delta * Math.PI / 180
  );
  createPolarizationVectors();
}

// Polarization mode switches
function switchToOMode() {
  conf.o_mode = true;
  conf.x_mode = false;
  updateFieldLines();
}

function switchToXMode() {
  conf.o_mode = false;
  conf.x_mode = true;
  updateFieldLines();
}

function toggleLockObserver() {
  if (conf.lock_observer) {
    controls.minPolarAngle = controls.maxPolarAngle = conf.obs_angle * Math.PI / 180;
  } else {
    controls.minPolarAngle = 0;
    controls.maxPolarAngle = Math.PI;
  }
}

// ============================================================================
// GUI CONTROLS
// ============================================================================

const gui = new GUI();
gui.add(conf, "obs_angle", 0.0, 90.0).listen().onChange(updatePolarizationSphere);
gui.add(conf, "pl_radius", 1.0, 10.0).listen().onChange(updatePolarizationSphere);
gui.add(conf, "dipole_angle", 0.0, 90.0).listen().onChange(updatePolarizationSphere);
gui.add(conf, "emission_h", 0.0, 1.0).listen().onChange(updatePolarizationSphere);
gui.add(conf, "translation_d", 0.0, 1.0).listen().onChange(updatePolarizationSphere);
gui.add(conf, "translation_delta", 0.0, 180.0).listen().onChange(updatePolarizationSphere);
gui.add(conf, "rotation_speed", 0.0, 1).listen();
gui.add(conf, "o_mode").listen().onChange(switchToOMode);
gui.add(conf, "x_mode").listen().onChange(switchToXMode);
gui.add(conf, "rotating").listen();
gui.add(conf, "lock_observer").listen().onChange(toggleLockObserver);

// Camera position setter
const guiFunctions = {
  set_observer: function () {
    camera.position.set(
      camera_distance * Math.sin(conf.obs_angle * Math.PI / 180),
      0,
      camera_distance * Math.cos(conf.obs_angle * Math.PI / 180)
    );
  }
};
gui.add(guiFunctions, "set_observer");

// ============================================================================
// ANIMATION LOOP
// ============================================================================

// Phase tracker for rotation
let phase = Math.PI / 2;

function animate3DScene() {
  requestAnimationFrame(animate3DScene);

  // Rotate field lines and polarization vectors if enabled
  if (conf.rotating) {
    const rotation_axis = new THREE.Vector3(0, 0, 1);
    field_lines.rotateOnWorldAxis(rotation_axis, -conf.rotation_speed);
    polarization_vectors.rotateOnWorldAxis(rotation_axis, -conf.rotation_speed);
    phase += conf.rotation_speed;
    phase = phase % (Math.PI * 4);
  }

  renderer.render(scene, camera);
  controls.update();
}

// Initialize field lines and start animation
createFieldLines(
  conf.obs_angle * Math.PI / 180,
  conf.pl_radius + conf.emission_h,
  40,
  100,
  conf.translation_d,
  conf.translation_delta * Math.PI / 180
);
createPolarizationVectors();
animate3DScene();

// ============================================================================
// 2D PLOTS - POLARIZATION ANGLE AND INTENSITY
// ============================================================================

// Canvas setup for 2D plots
const pa_canvas = document.getElementById("plot1");      // Polarization angle plot
const intensity_canvas = document.getElementById("plot2"); // Intensity plot
const devicePixelRatio = window.devicePixelRatio || 1;

// Plot dimensions
const numPoints_X = pa_canvas.width * 4;  // Number of points along x-axis
const numPoints_Y = pa_canvas.height;     // Number of points along y-axis

// ============================================================================
// COLOR DEFINITIONS
// ============================================================================

const colors = {
  red: new WebglPlotBundle.ColorRGBA(243/255, 91/255, 4/255, 1),     // Phase marker
  white: new WebglPlotBundle.ColorRGBA(1, 1, 1, 1),                   // Background PA
  green: new WebglPlotBundle.ColorRGBA(0, 1, 0, 1),                   // Observed PA
  blue: new WebglPlotBundle.ColorRGBA(26/255, 90/255, 217/255, 1),   // Intensity
  yellow: new WebglPlotBundle.ColorRGBA(1, 1, 0, 1)                   // (Unused - for reference)
};

// ============================================================================
// LINE DEFINITIONS
// ============================================================================

const line_thickness = 0.01;

// Polarization angle lines (split by intensity threshold)
const pa_line_observed = new WebglPlotBundle.WebglThickLine(colors.green, numPoints_X, line_thickness);
const pa_line_background = new WebglPlotBundle.WebglThickLine(colors.white, numPoints_X, line_thickness);

// Phase marker lines (vertical red line showing current rotation phase)
const phase_marker_pa = new WebglPlotBundle.WebglThickLine(colors.red, numPoints_Y, line_thickness);
const phase_marker_intensity = new WebglPlotBundle.WebglThickLine(colors.red, numPoints_Y, line_thickness);

// Intensity line
const intensity_line = new WebglPlotBundle.WebglThickLine(colors.blue, numPoints_X, line_thickness);

// ============================================================================
// WEBGL PLOT INITIALIZATION
// ============================================================================

// Create plot instances
const pa_plot = new WebglPlotBundle.WebglPlot(pa_canvas);           // PA plot
const intensity_plot = new WebglPlotBundle.WebglPlot(intensity_canvas); // Intensity plot

// Initialize X-axis spacing for all lines
pa_line_observed.lineSpaceX(-1, 2 / numPoints_X);
pa_line_background.lineSpaceX(-1, 2 / numPoints_X);
phase_marker_pa.lineSpaceX(-1, 2 / numPoints_Y);
intensity_line.lineSpaceX(-1, 2 / numPoints_X);
phase_marker_intensity.lineSpaceX(-1, 2 / numPoints_Y);

// Adjust intensity plot vertical offset (show range 0-2 instead of -1 to 1)
intensity_plot.gOffsetY = -0.5;

// Add lines to plots
pa_plot.addThickLine(pa_line_observed);
pa_plot.addThickLine(pa_line_background);
pa_plot.addThickLine(phase_marker_pa);

intensity_plot.addThickLine(intensity_line);
intensity_plot.addThickLine(phase_marker_intensity);

// ============================================================================
// PHYSICS CALCULATIONS
// ============================================================================

/**
 * Calculate Gaussian cone intensity profile
 * @param {number} theta_obs - Observer angle (radians)
 * @param {number} theta_cone - Cone opening angle (radians)
 * @param {number} dipole - Dipole angle (radians)
 * @param {number} phi - Rotation phase (radians)
 * @param {number} I0 - Peak intensity
 * @returns {number} Intensity at given position
 */
function calculateConeIntensity(theta_obs, theta_cone, dipole, phi, I0) {
  // Cone axis direction
  const ax = Math.sin(dipole) * Math.cos(phi);
  const ay = Math.sin(dipole) * Math.sin(phi);
  const az = Math.cos(dipole);

  // Observer position (on xz-plane)
  const r = 1.0;
  const rx = r * Math.sin(theta_obs);
  const ry = 0;
  const rz = r * Math.cos(theta_obs);

  // Perpendicular distance from cone axis
  const projection = rx * ax + rz * az;
  const r_perp_squared = Math.abs((rx * rx + rz * rz) - projection * projection);

  // Gaussian width parameter
  const sigma = 0.2 * 2 * Math.pow(Math.tan(theta_cone + 1e-8), 2);

  // Gaussian intensity profile
  return I0 * Math.exp(-r_perp_squared / sigma);
}

/**
 * Calculate polarization angle using RVM formula
 * @param {number} phi - Rotation phase (radians)
 * @param {number} dipole_angle - Dipole angle (degrees)
 * @param {number} obs_angle - Observer angle (degrees)
 * @param {number} emission_h - Emission height
 * @param {number} translation_d - Translation distance
 * @param {number} translation_delta - Translation angle (degrees)
 * @param {boolean} o_mode - True for O-mode, false for X-mode
 * @returns {number} Polarization angle (radians)
 */
function calculatePolarizationAngle(phi, dipole_angle, obs_angle, emission_h, translation_d, translation_delta, o_mode) {
  // Convert to radians
  const alpha = dipole_angle * Math.PI / 180;
  const xi = obs_angle * Math.PI / 180;
  const delta = translation_delta * Math.PI / 180;
  const beta = 90 / 180 * Math.PI;
  
  // RVM parameters
  const eta = emission_h;
  const epsilon = translation_d;

  // RVM equation components (Eq. 11 from RVM paper)
  const numerator = (1 + eta - epsilon * Math.cos(delta) * Math.cos(xi)) * Math.sin(alpha) * Math.sin(beta + phi)
                  + epsilon * Math.sin(delta) * (Math.cos(alpha) * Math.cos(xi) * Math.sin(phi) 
                  - Math.sin(alpha) * Math.sin(beta) * Math.sin(xi));
                  
  const denominator = (1 + eta) * (Math.cos(alpha) * Math.sin(xi) - Math.sin(alpha) * Math.cos(xi) * Math.cos(beta + phi))
                    + epsilon * (Math.sin(alpha) * Math.cos(delta) * Math.cos(beta + phi) 
                    - Math.cos(alpha) * Math.sin(delta) * Math.cos(phi));

  if (o_mode) {
    // O-mode (ordinary) polarization
    return Math.atan(numerator / denominator);
  } else {
    // X-mode (extraordinary) polarization
    let pa = Math.atan(numerator / denominator) + Math.PI / 2;
    if (pa > Math.PI / 2) {
      pa = pa - Math.PI;
    }
    return pa;
  }
}

// ============================================================================
// PLOT UPDATE FUNCTION
// ============================================================================

/**
 * Update all 2D plot lines based on current configuration
 */
function update2DPlots() {
  const noise_amplitude = 0.01;
  const pa_amplitude = 0.4;
  const camera_phi = Math.atan2(camera.position.y, camera.position.x) + Math.PI;
  
  // Intensity threshold for masking PA (1% of peak)
  const I0 = 1;
  const intensity_threshold = 0.01 * I0;
  
  // Pre-calculate all PA and intensity values
  const intensity_array = new Array(numPoints_X);
  const pa_array = new Array(numPoints_X);
  
  for (let i = 0; i < numPoints_X; i++) {
    const x = -1 + i * 2 / numPoints_X;
    const phi = 2.0 * Math.PI * (x + 1);
    
    // Calculate polarization angle
    const pa = calculatePolarizationAngle(
      phi,
      conf.dipole_angle,
      conf.obs_angle,
      conf.emission_h,
      conf.translation_d,
      conf.translation_delta,
      conf.o_mode
    );
    pa_array[i] = pa * pa_amplitude + (Math.random() - 0.5) * noise_amplitude;
    
    // Calculate intensity
    const theta = conf.obs_angle * Math.PI / 180;
    const theta_cone = Math.PI / 18;  // Fixed cone angle
    const dipole = conf.dipole_angle * Math.PI / 180;
    intensity_array[i] = calculateConeIntensity(theta, theta_cone, dipole, phi, I0);
  }
  
  // Set PA lines based on intensity threshold
  for (let i = 0; i < numPoints_X; i++) {
    if (intensity_array[i] > intensity_threshold) {
      // High intensity region: show observed PA (green), hide background
      pa_line_observed.setY(i, pa_array[i]);
      pa_line_background.setY(i, NaN);
    } else {
      // Low intensity region: show background PA (white), hide observed
      pa_line_observed.setY(i, NaN);
      pa_line_background.setY(i, pa_array[i]);
    }
  }

  // Update phase marker on PA plot (vertical red line)
  for (let i = 0; i < phase_marker_pa.numPoints; i++) {
    const angle = (camera_phi + phase) % (Math.PI * 4);
    const x = angle / (Math.PI * 4) * 2 - 1;  // Scale to [-1, 1] range
    phase_marker_pa.setX(i, x);
    phase_marker_pa.setY(i, -1 + i * 2 / numPoints_Y);
  }

  // Update intensity line
  for (let i = 0; i < intensity_line.numPoints; i++) {
    intensity_line.setY(i, intensity_array[i] + (Math.random() - 0.5) * noise_amplitude);
  }

  // Update phase marker on intensity plot
  for (let i = 0; i < phase_marker_intensity.numPoints; i++) {
    const angle = (camera_phi + phase) % (Math.PI * 4);
    const x = angle / (Math.PI * 4) * 2 - 1;
    phase_marker_intensity.setX(i, x);
    phase_marker_intensity.setY(i, -1 + i * 2 / numPoints_Y + 0.5);
  }
}

// ============================================================================
// 2D PLOT ANIMATION LOOP
// ============================================================================

function animate2DPlots() {
  update2DPlots();
  pa_plot.update();
  intensity_plot.update();
  requestAnimationFrame(animate2DPlots);
}

animate2DPlots();

// ============================================================================
// WEBGL PLOT BUNDLE (External Library)
// ============================================================================
// Note: The WebglPlotBundle is included inline here for convenience.
// In production, this should be loaded as a separate module.

(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
  typeof define === 'function' && define.amd ? define(['exports'], factory) :
  (global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.WebglPlotBundle = {}));
})(this, (function (exports) { 'use strict';

  class ColorRGBA {
      r; g; b; a;
      constructor(r, g, b, a) {
          this.r = r; this.g = g; this.b = b; this.a = a;
      }
  }

  class WebglBase {
      intensity; visible; numPoints; xy; color;
      scaleX; scaleY; offsetX; offsetY; loop;
      webglNumPoints; _vbuffer; _coord;
      
      constructor() {
          this.scaleX = 1; this.scaleY = 1;
          this.offsetX = 0; this.offsetY = 0;
          this.loop = false; this._vbuffer = 0; this._coord = 0;
          this.visible = true; this.intensity = 1;
          this.xy = new Float32Array([]);
          this.numPoints = 0;
          this.color = new ColorRGBA(0, 0, 0, 1);
          this.webglNumPoints = 0;
      }
  }

  class WebglLine extends WebglBase {
      currentIndex = 0;
      constructor(c, numPoints) {
          super();
          this.webglNumPoints = numPoints;
          this.numPoints = numPoints;
          this.color = c;
          this.xy = new Float32Array(2 * this.webglNumPoints);
      }
      setX(index, x) { this.xy[index * 2] = x; }
      setY(index, y) { this.xy[index * 2 + 1] = y; }
      getX(index) { return this.xy[index * 2]; }
      getY(index) { return this.xy[index * 2 + 1]; }
      lineSpaceX(start, stepSize) {
          for (let i = 0; i < this.numPoints; i++) {
              this.setX(i, start + stepSize * i);
          }
      }
      arrangeX() { this.lineSpaceX(-1, 2 / this.numPoints); }
      constY(c) {
          for (let i = 0; i < this.numPoints; i++) { this.setY(i, c); }
      }
      shiftAdd(data) {
          const shiftSize = data.length;
          for (let i = 0; i < this.numPoints - shiftSize; i++) {
              this.setY(i, this.getY(i + shiftSize));
          }
          for (let i = 0; i < shiftSize; i++) {
              this.setY(i + this.numPoints - shiftSize, data[i]);
          }
      }
      addArrayY(yArray) {
          if (this.currentIndex + yArray.length <= this.numPoints) {
              for (let i = 0; i < yArray.length; i++) {
                  this.setY(this.currentIndex, yArray[i]);
                  this.currentIndex++;
              }
          }
      }
      replaceArrayY(yArray) {
          if (yArray.length == this.numPoints) {
              for (let i = 0; i < this.numPoints; i++) {
                  this.setY(i, yArray[i]);
              }
          }
      }
  }

  class WebglStep extends WebglBase {
      constructor(c, num) {
          super();
          this.webglNumPoints = num * 2;
          this.numPoints = num;
          this.color = c;
          this.xy = new Float32Array(2 * this.webglNumPoints);
      }
      setY(index, y) {
          this.xy[index * 4 + 1] = y;
          this.xy[index * 4 + 3] = y;
      }
      getX(index) { return this.xy[index * 4]; }
      getY(index) { return this.xy[index * 4 + 1]; }
      lineSpaceX(start, stepsize) {
          for (let i = 0; i < this.numPoints; i++) {
              this.xy[i * 4] = start + i * stepsize;
              this.xy[i * 4 + 2] = start + (i * stepsize + stepsize);
          }
      }
      constY(c) {
          for (let i = 0; i < this.numPoints; i++) { this.setY(i, c); }
      }
      shiftAdd(data) {
          const shiftSize = data.length;
          for (let i = 0; i < this.numPoints - shiftSize; i++) {
              this.setY(i, this.getY(i + shiftSize));
          }
          for (let i = 0; i < shiftSize; i++) {
              this.setY(i + this.numPoints - shiftSize, data[i]);
          }
      }
  }

  class WebglPolar extends WebglBase {
      numPoints; xy; color; intenisty; visible; offsetTheta;
      constructor(c, numPoints) {
          super();
          this.webglNumPoints = numPoints; this.numPoints = numPoints;
          this.color = c; this.intenisty = 1;
          this.xy = new Float32Array(2 * this.webglNumPoints);
          this._vbuffer = 0; this._coord = 0;
          this.visible = true; this.offsetTheta = 0;
      }
      setRtheta(index, theta, r) {
          const x = r * Math.cos((2 * Math.PI * (theta + this.offsetTheta)) / 360);
          const y = r * Math.sin((2 * Math.PI * (theta + this.offsetTheta)) / 360);
          this.setX(index, x);
          this.setY(index, y);
      }
      getTheta(index) { return 0; }
      getR(index) {
          return Math.sqrt(Math.pow(this.getX(index), 2) + Math.pow(this.getY(index), 2));
      }
      setX(index, x) { this.xy[index * 2] = x; }
      setY(index, y) { this.xy[index * 2 + 1] = y; }
      getX(index) { return this.xy[index * 2]; }
      getY(index) { return this.xy[index * 2 + 1]; }
  }

  class WebglSquare extends WebglBase {
      constructor(c) {
          super();
          this.webglNumPoints = 4; this.numPoints = 4;
          this.color = c;
          this.xy = new Float32Array(2 * this.webglNumPoints);
      }
      setSquare(x1, y1, x2, y2) {
          this.xy = new Float32Array([x1, y1, x1, y2, x2, y1, x2, y2]);
      }
  }

  // Vector math utilities
  const scaleAndAdd = (a, b, scale) => {
      return { x: a.x + b.x * scale, y: a.y + b.y * scale };
  };
  const normal = (dir) => {
      return { x: -dir.y, y: dir.x };
  };
  const direction = (a, b) => {
      let out = subtract(a, b);
      return normalize(out);
  };
  const add = (a, b) => {
      return { x: a.x + b.x, y: a.y + b.y };
  };
  const dot = (a, b) => {
      return a.x * b.x + a.y * b.y;
  };
  const normalize = (a) => {
      const out = { x: 0, y: 0 };
      let len = a.x * a.x + a.y * a.y;
      if (len > 0) {
          len = 1 / Math.sqrt(len);
          out.x = a.x * len;
          out.y = a.y * len;
      }
      return out;
  };
  const set = (x, y) => {
      return { x: x, y: y };
  };
  const subtract = (a, b) => {
      return { x: a.x - b.x, y: a.y - b.y };
  };

  const PolyLine = (lineXY) => {
      let curNormal, lineA = { x: 0, y: 0 }, lineB = { x: 0, y: 0 };
      const out = [];
      const addNext = (normal, length) => {
          out.push({ vec2: normal, miterLength: length });
      };
      const getXY = (index) => {
          return { x: lineXY[index * 2], y: lineXY[index * 2 + 1] };
      };
      lineA = direction(getXY(1), getXY(0));
      curNormal = normal(lineA);
      addNext(curNormal, 1);
      const numPoints = lineXY.length / 2;
      for (let i = 1; i < numPoints - 1; i++) {
          const last = getXY(i - 1);
          const cur = getXY(i);
          const next = getXY(i + 1);
          lineA = direction(cur, last);
          curNormal = normal(lineA);
          lineB = direction(next, cur);
          const miter = computeMiter(lineA, lineB);
          const miterLen = computeMiterLen(lineA, miter, 1);
          addNext(miter, miterLen);
      }
      lineA = direction(getXY(numPoints - 1), getXY(numPoints - 2));
      curNormal = normal(lineA);
      addNext(curNormal, 1);
      return out;
  };
  const computeMiter = (lineA, lineB) => {
      let tangent = add(lineA, lineB);
      tangent = normalize(tangent);
      return set(-tangent.y, tangent.x);
  };
  const computeMiterLen = (lineA, miter, halfThick) => {
      const tmp = set(-lineA.y, lineA.x);
      return halfThick / dot(miter, tmp);
  };

  class WebglThickLine extends WebglBase {
      currentIndex = 0;
      _linePoints; _thicknessRequested = 0; _actualThickness = 0;
      
      constructor(c, numPoints, thickness) {
          super();
          this.webglNumPoints = numPoints * 2;
          this.numPoints = numPoints;
          this.color = c;
          this._thicknessRequested = thickness;
          this._linePoints = new Float32Array(numPoints * 2);
          this.xy = new Float32Array(2 * this.webglNumPoints);
      }
      convertToTriPoints() {
          const halfThick = this._actualThickness / 2;
          const normals = PolyLine(this._linePoints);
          for (let i = 0; i < this.numPoints; i++) {
              const x = this._linePoints[2 * i];
              const y = this._linePoints[2 * i + 1];
              const point = { x: x, y: y };
              const top = scaleAndAdd(point, normals[i].vec2, normals[i].miterLength * halfThick);
              const bot = scaleAndAdd(point, normals[i].vec2, -normals[i].miterLength * halfThick);
              this.xy[i * 4] = top.x;
              this.xy[i * 4 + 1] = top.y;
              this.xy[i * 4 + 2] = bot.x;
              this.xy[i * 4 + 3] = bot.y;
          }
      }
      setX(index, x) { this._linePoints[index * 2] = x; }
      setY(index, y) { this._linePoints[index * 2 + 1] = y; }
      lineSpaceX(start, stepSize) {
          for (let i = 0; i < this.numPoints; i++) {
              this.setX(i, start + stepSize * i);
          }
      }
      setThickness(thickness) { this._thicknessRequested = thickness; }
      getThickness() { return this._thicknessRequested; }
      setActualThickness(thickness) { this._actualThickness = thickness; }
  }

  class WebglPlot {
      webgl; gScaleX; gScaleY; gXYratio; gOffsetX; gOffsetY;
      gLog10X; gLog10Y; _linesData; _linesAux; _thickLines; _surfaces;
      _progLine; debug = false;
      
      get linesData() { return this._linesData; }
      get linesAux() { return this._linesAux; }
      get thickLines() { return this._thickLines; }
      get surfaces() { return this._surfaces; }
      
      constructor(canvas, options) {
          if (options == undefined) {
              this.webgl = canvas.getContext("webgl", {
                  antialias: true,
                  transparent: false,
              });
          }
          else {
              this.webgl = canvas.getContext("webgl", {
                  antialias: options.antialias,
                  transparent: options.transparent,
                  desynchronized: options.deSync,
                  powerPerformance: options.powerPerformance,
                  preserveDrawing: options.preserveDrawing,
              });
              this.debug = options.debug == undefined ? false : options.debug;
          }
          this.log("canvas type is: " + canvas.constructor.name);
          this.log(`[webgl-plot]:width=${canvas.width}, height=${canvas.height}`);
          this._linesData = []; this._linesAux = [];
          this._thickLines = []; this._surfaces = [];
          this.gScaleX = 1; this.gScaleY = 1; this.gXYratio = 1;
          this.gOffsetX = 0; this.gOffsetY = 0;
          this.gLog10X = false; this.gLog10Y = false;
          this.webgl.clear(this.webgl.COLOR_BUFFER_BIT);
          this.webgl.viewport(0, 0, canvas.width, canvas.height);
          this._progLine = this.webgl.createProgram();
          this.initThinLineProgram();
          this.webgl.enable(this.webgl.BLEND);
          this.webgl.blendFunc(this.webgl.SRC_ALPHA, this.webgl.ONE_MINUS_SRC_ALPHA);
      }
      _drawLines(lines) {
          const webgl = this.webgl;
          lines.forEach((line) => {
              if (line.visible) {
                  webgl.useProgram(this._progLine);
                  const uscale = webgl.getUniformLocation(this._progLine, "uscale");
                  webgl.uniformMatrix2fv(uscale, false, new Float32Array([
                      line.scaleX * this.gScaleX * (this.gLog10X ? 1 / Math.log(10) : 1), 0, 0,
                      line.scaleY * this.gScaleY * this.gXYratio * (this.gLog10Y ? 1 / Math.log(10) : 1),
                  ]));
                  const uoffset = webgl.getUniformLocation(this._progLine, "uoffset");
                  webgl.uniform2fv(uoffset, new Float32Array([line.offsetX + this.gOffsetX, line.offsetY + this.gOffsetY]));
                  const isLog = webgl.getUniformLocation(this._progLine, "is_log");
                  webgl.uniform2iv(isLog, new Int32Array([this.gLog10X ? 1 : 0, this.gLog10Y ? 1 : 0]));
                  const uColor = webgl.getUniformLocation(this._progLine, "uColor");
                  webgl.uniform4fv(uColor, [line.color.r, line.color.g, line.color.b, line.color.a]);
                  webgl.bufferData(webgl.ARRAY_BUFFER, line.xy, webgl.STREAM_DRAW);
                  webgl.drawArrays(line.loop ? webgl.LINE_LOOP : webgl.LINE_STRIP, 0, line.webglNumPoints);
              }
          });
      }
      _drawSurfaces(squares) {
          const webgl = this.webgl;
          squares.forEach((square) => {
              if (square.visible) {
                  webgl.useProgram(this._progLine);
                  const uscale = webgl.getUniformLocation(this._progLine, "uscale");
                  webgl.uniformMatrix2fv(uscale, false, new Float32Array([
                      square.scaleX * this.gScaleX * (this.gLog10X ? 1 / Math.log(10) : 1), 0, 0,
                      square.scaleY * this.gScaleY * this.gXYratio * (this.gLog10Y ? 1 / Math.log(10) : 1),
                  ]));
                  const uoffset = webgl.getUniformLocation(this._progLine, "uoffset");
                  webgl.uniform2fv(uoffset, new Float32Array([square.offsetX + this.gOffsetX, square.offsetY + this.gOffsetY]));
                  const isLog = webgl.getUniformLocation(this._progLine, "is_log");
                  webgl.uniform2iv(isLog, new Int32Array([this.gLog10X ? 1 : 0, this.gLog10Y ? 1 : 0]));
                  const uColor = webgl.getUniformLocation(this._progLine, "uColor");
                  webgl.uniform4fv(uColor, [square.color.r, square.color.g, square.color.b, square.color.a]);
                  webgl.bufferData(webgl.ARRAY_BUFFER, square.xy, webgl.STREAM_DRAW);
                  webgl.drawArrays(webgl.TRIANGLE_STRIP, 0, square.webglNumPoints);
              }
          });
      }
      _drawTriangles(thickLine) {
          const webgl = this.webgl;
          webgl.bufferData(webgl.ARRAY_BUFFER, thickLine.xy, webgl.STREAM_DRAW);
          webgl.useProgram(this._progLine);
          const uscale = webgl.getUniformLocation(this._progLine, "uscale");
          webgl.uniformMatrix2fv(uscale, false, new Float32Array([
              thickLine.scaleX * this.gScaleX * (this.gLog10X ? 1 / Math.log(10) : 1), 0, 0,
              thickLine.scaleY * this.gScaleY * this.gXYratio * (this.gLog10Y ? 1 / Math.log(10) : 1),
          ]));
          const uoffset = webgl.getUniformLocation(this._progLine, "uoffset");
          webgl.uniform2fv(uoffset, new Float32Array([thickLine.offsetX + this.gOffsetX, thickLine.offsetY + this.gOffsetY]));
          const isLog = webgl.getUniformLocation(this._progLine, "is_log");
          webgl.uniform2iv(isLog, new Int32Array([0, 0]));
          const uColor = webgl.getUniformLocation(this._progLine, "uColor");
          webgl.uniform4fv(uColor, [
              thickLine.color.r, thickLine.color.g,
              thickLine.color.b, thickLine.color.a,
          ]);
          webgl.drawArrays(webgl.TRIANGLE_STRIP, 0, thickLine.xy.length / 2);
      }
      _drawThickLines() {
          this._thickLines.forEach((thickLine) => {
              if (thickLine.visible) {
                  const calibFactor = Math.min(this.gScaleX, this.gScaleY);
                  thickLine.setActualThickness(thickLine.getThickness() / calibFactor);
                  thickLine.convertToTriPoints();
                  this._drawTriangles(thickLine);
              }
          });
      }
      update() {
          this.clear();
          this.draw();
      }
      draw() {
          this._drawLines(this.linesData);
          this._drawLines(this.linesAux);
          this._drawThickLines();
          this._drawSurfaces(this.surfaces);
      }
      clear() {
          this.webgl.clear(this.webgl.COLOR_BUFFER_BIT);
      }
      _addLine(line) {
          line._vbuffer = this.webgl.createBuffer();
          this.webgl.bindBuffer(this.webgl.ARRAY_BUFFER, line._vbuffer);
          this.webgl.bufferData(this.webgl.ARRAY_BUFFER, line.xy, this.webgl.STREAM_DRAW);
          line._coord = this.webgl.getAttribLocation(this._progLine, "coordinates");
          this.webgl.vertexAttribPointer(line._coord, 2, this.webgl.FLOAT, false, 0, 0);
          this.webgl.enableVertexAttribArray(line._coord);
      }
      addDataLine(line) {
          this._addLine(line);
          this.linesData.push(line);
      }
      addLine = this.addDataLine;
      addAuxLine(line) {
          this._addLine(line);
          this.linesAux.push(line);
      }
      addThickLine(thickLine) {
          this._addLine(thickLine);
          this._thickLines.push(thickLine);
      }
      addSurface(surface) {
          this._addLine(surface);
          this.surfaces.push(surface);
      }
      initThinLineProgram() {
          const vertCode = `
    attribute vec2 coordinates;
    uniform mat2 uscale;
    uniform vec2 uoffset;
    uniform ivec2 is_log;
    void main(void) {
       float x = (is_log[0]==1) ? log(coordinates.x) : coordinates.x;
       float y = (is_log[1]==1) ? log(coordinates.y) : coordinates.y;
       vec2 line = vec2(x, y);
       gl_Position = vec4(uscale*line + uoffset, 0.0, 1.0);
    }`;
          const vertShader = this.webgl.createShader(this.webgl.VERTEX_SHADER);
          this.webgl.shaderSource(vertShader, vertCode);
          this.webgl.compileShader(vertShader);
          const fragCode = `
       precision mediump float;
       uniform highp vec4 uColor;
       void main(void) {
          gl_FragColor =  uColor;
       }`;
          const fragShader = this.webgl.createShader(this.webgl.FRAGMENT_SHADER);
          this.webgl.shaderSource(fragShader, fragCode);
          this.webgl.compileShader(fragShader);
          this._progLine = this.webgl.createProgram();
          this.webgl.attachShader(this._progLine, vertShader);
          this.webgl.attachShader(this._progLine, fragShader);
          this.webgl.linkProgram(this._progLine);
      }
      popDataLine() { this.linesData.pop(); }
      removeAllLines() {
          this._linesData = []; this._linesAux = [];
          this._thickLines = []; this._surfaces = [];
      }
      removeDataLines() { this._linesData = []; }
      removeAuxLines() { this._linesAux = []; }
      viewport(a, b, c, d) {
          this.webgl.viewport(a, b, c, d);
      }
      log(str) {
          if (this.debug) {
              console.log("[webgl-plot]:" + str);
          }
      }
  }

  exports.ColorRGBA = ColorRGBA;
  exports.WebglLine = WebglLine;
  exports.WebglPlot = WebglPlot;
  exports.WebglPolar = WebglPolar;
  exports.WebglSquare = WebglSquare;
  exports.WebglStep = WebglStep;
  exports.WebglThickLine = WebglThickLine;

}));

// ============================================================================
// UNUSED CODE (Commented out for reference)
// ============================================================================

/*
// Chart.js scatter plot (unused)
const xyValues = [
    {x:0, y:-1},
    {x:0, y:1}
];

const color1 = new WebglPlotBundle.ColorRGBA(0,0.1,0.1,0.5);

new Chart("chart", {
    type: "scatter",
    data: {},
    options: {
        scales: {
            xAxes: [{ticks: {min: 0, max:1},display:false}],
            yAxes: [{ticks: {min: -1, max:1}, display:true,gridLines: {
                display: true         
            }}],
        },
        responsive: false,
    }
});
*/
