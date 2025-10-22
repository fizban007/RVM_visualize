/**
 * ===================================================================
 * ROTATING VECTOR MODEL (RVM) VISUALIZATION
 * ===================================================================
 * 
 * This script visualizes the rotating vector model of pulsar emission,
 * showing magnetic field lines, polarization vectors, and real-time
 * polarization angle calculations.
 * 
 * Main Components:
 * - 3D visualization using THREE.js
 * - Real-time plotting using WebGL
 * - Interactive controls using dat.GUI
 * - Physics calculations for RVM
 */

// ===================================================================
// IMPORTS AND DEPENDENCIES
// ===================================================================
import * as THREE from '../node_modules/three/build/three.module.js';
import { GUI } from '../node_modules/dat.gui/build/dat.gui.module.js';
import { OrbitControls } from '../node_modules/three/examples/jsm/controls/OrbitControls.js';
import { Line2 } from '../node_modules/three/examples/jsm/lines/Line2.js';
import { LineMaterial } from '../node_modules/three/examples/jsm/lines/LineMaterial.js';
import { LineGeometry } from '../node_modules/three/examples/jsm/lines/LineGeometry.js';

// ===================================================================
// GLOBAL CONFIGURATION
// ===================================================================

/**
 * Configuration object containing all physical and visualization parameters
 */
class PulsarConfig {
  constructor() {
    // Physical parameters
    this.dipole_angle = 30.0;        // Angle between magnetic and rotation axes (degrees)
    this.obs_angle = 90.0;           // Observer's viewing angle from rotation axis (degrees)
    this.pl_radius = 10.0;           // Radius of polarization limiting sphere
    this.emission_h = 0.01;          // Height of radio emission region
    this.translation_d = 0.41;       // Offset distance of emission beam
    this.translation_delta = 0.11;   // Angular aperture of emission beam (radians)
    
    // Animation parameters
    this.rotation_speed = 0.01;      // Speed of pulsar rotation
    
    // Display modes
    this.x_mode = false;             // Show X-mode (extraordinary) polarization
    this.o_mode = true;              // Show O-mode (ordinary) polarization
    this.rotating = false;           // Enable/disable rotation animation
    this.lock_observer = true;       // Lock observer viewpoint
  }
}

// Global variables
const config = new PulsarConfig();
const quaternions = {
  dipole: new THREE.Quaternion(),
  dipole_inv: new THREE.Quaternion()
};

// Canvas and rendering setup
const canvas = document.getElementById("vis");
const width = window.innerWidth;
const height = window.innerHeight;

// ===================================================================
// RENDERER SETUP
// ===================================================================

const renderer = new THREE.WebGLRenderer({
  canvas: canvas,
  alpha: true,
  antialias: true,
  preserveDrawingBuffer: true
});
renderer.setSize(width, height);

// ===================================================================
// SCENE SETUP
// ===================================================================

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x000000);

// Camera setup
const camera = new THREE.PerspectiveCamera(30, width / height, 1, 1000);
const camera_distance = 160;
camera.position.z = camera_distance * Math.cos(80 * Math.PI / 180);
camera.position.x = camera_distance * Math.sin(80 * Math.PI / 180);
camera.up.set(0, 0, 1);
camera.lookAt(new THREE.Vector3(0, 0, 0));

// Controls setup
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.05;
controls.enableKeys = false;
controls.minPolarAngle = controls.maxPolarAngle = config.obs_angle * Math.PI / 180;

// ===================================================================
// LIGHTING SETUP
// ===================================================================

function setupLighting() {
  // Directional light
  const directionalLight = new THREE.DirectionalLight(0xffffff);
  directionalLight.position.set(107, 107, 107);
  scene.add(directionalLight);

  // Ambient light
  const ambientLight = new THREE.AmbientLight(0x404040);
  scene.add(ambientLight);
}

// ===================================================================
// 3D OBJECTS CREATION
// ===================================================================

/**
 * Creates the neutron star sphere
 */
function createNeutronStar() {
  const geometry = new THREE.SphereGeometry(1.0, 64, 64);
  const material = new THREE.MeshPhongMaterial({ color: 0xaaaaaa });
  const mesh = new THREE.Mesh(geometry, material);
  scene.add(mesh);
  return mesh;
}

/**
 * Creates the polarization limiting sphere
 */
function createPolarizationSphere() {
  const geometry = new THREE.SphereGeometry(1.0, 64, 64);
  const material = new THREE.MeshPhongMaterial({
    color: 0xffffff,
    transparent: true,
    depthWrite: false,
    blending: THREE.NormalBlending,
    opacity: 0.3
  });
  const mesh = new THREE.Mesh(geometry, material);
  mesh.scale.setScalar(config.pl_radius);
  scene.add(mesh);
  return mesh;
}

/**
 * Creates the rotation axis line
 */
function createRotationAxis() {
  const geometry = new THREE.BufferGeometry();
  const vertices = new Float32Array([0, 0, -300, 0, 0, 300]);
  geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
  
  const material = new THREE.LineBasicMaterial({
    color: 0x8080aa,
    linewidth: 2.5
  });
  
  const line = new THREE.Line(geometry, material);
  scene.add(line);
  return line;
}

/**
 * Creates the magnetic dipole beams (cones)
 */
function createMagneticBeams() {
  const loader = new THREE.TextureLoader();
  const texture = loader.load('plasma.jpg');
  
  const material = new THREE.MeshPhongMaterial({
    map: texture,
    transparent: true,
    opacity: 0.7,
    blending: THREE.NormalBlending,
    depthWrite: false,
    side: THREE.DoubleSide,
    alphaTest: 0.1
  });

  const fieldLines = new THREE.Group();

  // Create two opposing cones for dipole field
  const geometry1 = new THREE.ConeGeometry(5, 100, 64, 32, true);
  const cone1 = new THREE.Mesh(geometry1, material);
  cone1.name = "line";
  cone1.rotateX(-(config.dipole_angle / 180 * Math.PI + Math.PI / 2));
  geometry1.translate(0, -50, 0);
  fieldLines.add(cone1);

  const geometry2 = new THREE.ConeGeometry(5, 100, 64, 32, true);
  const cone2 = new THREE.Mesh(geometry2, material);
  cone2.name = "line";
  cone2.rotateX(-(config.dipole_angle / 180 * Math.PI - Math.PI / 2));
  geometry2.translate(0, -50, 0);
  fieldLines.add(cone2);

  scene.add(fieldLines);
  return fieldLines;
}

// ===================================================================
// PHYSICS CALCULATIONS
// ===================================================================

/**
 * Updates quaternions based on current dipole angle
 */
function updateQuaternions() {
  quaternions.dipole.setFromAxisAngle(
    new THREE.Vector3(1, 0, 0), 
    Math.PI * config.dipole_angle / 180.0
  );
  quaternions.dipole_inv.copy(quaternions.dipole);
  quaternions.dipole_inv.invert();
}

/**
 * Calculates field lines based on observer angle and other parameters
 */
function calculateFieldLines(th_obs, r_pl, n_lines, n_segments, d, delta) {
  const intersectPoints = [];
  
  for (let i = 0; i < n_lines; i++) {
    const phi = i * 2.0 * Math.PI / n_lines;
    
    // Calculate offset
    const offset = new THREE.Vector3(
      d * Math.sin(delta) * Math.cos(phi),
      d * Math.sin(delta) * Math.sin(phi),
      d * Math.cos(delta)
    );

    // Calculate intersection point in lab frame
    const p_intersect_lab = new THREE.Vector3(
      r_pl * Math.sin(th_obs) * Math.cos(phi),
      r_pl * Math.sin(th_obs) * Math.sin(phi),
      r_pl * Math.cos(th_obs)
    );
    
    intersectPoints.push(p_intersect_lab.clone());
    
    // Transform to magnetic dipole frame
    const p_intersect = p_intersect_lab.clone().sub(offset);
    p_intersect.applyQuaternion(quaternions.dipole);
    
    // Additional field line calculations would go here...
  }
  
  return intersectPoints;
}

// ===================================================================
// WEBGL PLOTTING SETUP
// ===================================================================

/**
 * Sets up the WebGL plotting canvas and lines
 */
function setupWebGLPlotting() {
  const plotCanvas = document.getElementById("plot");
  
  // Canvas setup
  const devicePixelRatio = window.devicePixelRatio || 1;
  plotCanvas.width = plotCanvas.clientWidth * devicePixelRatio;
  plotCanvas.height = plotCanvas.clientHeight * devicePixelRatio * 1.4;
  
  const numX = plotCanvas.width;
  const numY = plotCanvas.height;
  
  // Color definitions
  const colors = {
    red: new WebglPlotBundle.ColorRGBA(243/255, 91/255, 4/255, 1),
    white: new WebglPlotBundle.ColorRGBA(1, 1, 1, 1),
    green: new WebglPlotBundle.ColorRGBA(0, 1, 0, 1)
  };
  
  const thickness = 0.01;
  
  // Create plot lines
  const lines = {
    polarization: new WebglPlotBundle.WebglThickLine(colors.white, numX, thickness),
    phase: new WebglPlotBundle.WebglThickLine(colors.red, numY, thickness),
    reference: new WebglPlotBundle.WebglThickLine(colors.green, numX, thickness)
  };
  
  // WebGL plot setup
  const wglp = new WebglPlotBundle.WebglPlot(plotCanvas);
  
  // Configure line spacing
  lines.polarization.lineSpaceX(-1, 2 / numX);
  lines.phase.lineSpaceX(-1, 2 / numX);
  lines.reference.lineSpaceX(-1, 2 / numX);
  
  // Add lines to plot
  wglp.addThickLine(lines.polarization);
  wglp.addThickLine(lines.phase);
  wglp.addThickLine(lines.reference);
  
  return { wglp, lines, numX, numY };
}

// ===================================================================
// ANIMATION AND UPDATE FUNCTIONS
// ===================================================================

let phase = 0;
const plotting = setupWebGLPlotting();

/**
 * Updates the WebGL plot data
 */
function updatePlotData() {
  const { lines, numX } = plotting;
  const freq = 0.002;
  const amp = 0.4;
  const noise = 0.00;
  
  // Calculate current phase
  const cam_phi = Math.atan2(camera.position.y, camera.position.x) + Math.PI;
  const angle = (cam_phi + phase) % (Math.PI * 4);
  const x = angle / (Math.PI * 4) * 2 - 1;
  
  // Update polarization line (line 1)
  for (let i = 0; i < lines.polarization.numPoints; i++) {
    const rad_dipole = config.dipole_angle / 180 * Math.PI;
    const rad_obs = config.obs_angle / 180 * Math.PI;
    const x = -1 + i * 2 / numX;
    const phi = 2.0 * Math.PI * (x + 1);
    
    // Physical parameters
    const eta = config.emission_h;
    const eps = config.translation_d;
    const xi = rad_obs;
    const delta = config.translation_delta * Math.PI / 180;
    const beta = 90 / 180 * Math.PI;
    
    // RVM equations
    const eq11a = (1 + eta - eps * Math.cos(delta) * Math.cos(xi)) * Math.sin(rad_dipole) * Math.sin(beta + phi)
                  + eps * Math.sin(delta) * (Math.cos(rad_dipole) * Math.cos(xi) * Math.sin(phi) 
                  - Math.sin(rad_dipole) * Math.sin(beta) * Math.sin(xi));
                  
    const eq11b = (1 + eta) * (Math.cos(rad_dipole) * Math.sin(xi) 
                  - Math.sin(rad_dipole) * Math.cos(xi) * Math.cos(beta + phi))
                  + eps * (Math.sin(rad_dipole) * Math.cos(delta) * Math.cos(beta + phi) 
                  - Math.cos(rad_dipole) * Math.sin(delta) * Math.cos(phi));
    
    // Calculate polarization angle
    let ySin;
    if (config.o_mode) {
      ySin = Math.atan(eq11a / eq11b);
    } else {
      ySin = Math.atan(eq11a / eq11b) + Math.PI / 2;
      if (ySin > Math.PI / 2) {
        ySin = ySin - Math.PI;
      }
    }
    
    const yNoise = Math.random() - 0.5;
    lines.polarization.setY(i, ySin * amp + yNoise * noise);
  }
  
  // Update phase marker line (line 2)
  for (let i = 0; i < lines.phase.numPoints; i++) {
    lines.phase.setX(i, x);
    lines.phase.setY(i, -1 + i * 2 / numX);
  }
}

/**
 * Main animation loop
 */
function animate() {
  // Update controls
  controls.update();
  
  // Update plot data
  updatePlotData();
  plotting.wglp.update();
  
  // Render 3D scene
  renderer.render(scene, camera);
  
  // Continue animation
  requestAnimationFrame(animate);
}

/**
 * Animation frame for WebGL plotting
 */
function plotAnimationFrame() {
  updatePlotData();
  plotting.wglp.update();
  requestAnimationFrame(plotAnimationFrame);
}

// ===================================================================
// GUI SETUP
// ===================================================================

/**
 * Sets up the dat.GUI interface
 */
function setupGUI() {
  const gui = new GUI();
  
  // Physical parameters folder
  const physicsFolder = gui.addFolder('Physics Parameters');
  physicsFolder.add(config, 'dipole_angle', 0, 90).name('Dipole Angle (°)').onChange(updateQuaternions);
  physicsFolder.add(config, 'obs_angle', 0, 180).name('Observer Angle (°)');
  physicsFolder.add(config, 'pl_radius', 1, 20).name('Polarization Radius');
  physicsFolder.add(config, 'emission_h', 0, 1).name('Emission Height');
  physicsFolder.add(config, 'translation_d', 0, 1).name('Beam Offset');
  physicsFolder.add(config, 'translation_delta', 0, 1).name('Beam Aperture');
  
  // Animation folder
  const animationFolder = gui.addFolder('Animation');
  animationFolder.add(config, 'rotation_speed', 0, 0.1).name('Rotation Speed');
  animationFolder.add(config, 'rotating').name('Enable Rotation');
  
  // Display modes folder
  const displayFolder = gui.addFolder('Display Modes');
  displayFolder.add(config, 'x_mode').name('X-mode Polarization');
  displayFolder.add(config, 'o_mode').name('O-mode Polarization');
  displayFolder.add(config, 'lock_observer').name('Lock Observer');
  
  // Open important folders by default
  physicsFolder.open();
  displayFolder.open();
  
  return gui;
}

// ===================================================================
// INITIALIZATION
// ===================================================================

/**
 * Initialize the entire application
 */
function initialize() {
  console.log('Initializing RVM Visualization...');
  
  // Setup scene components
  setupLighting();
  createNeutronStar();
  createPolarizationSphere();
  createRotationAxis();
  createMagneticBeams();
  
  // Initialize physics
  updateQuaternions();
  
  // Setup GUI
  setupGUI();
  
  // Start animations
  animate();
  plotAnimationFrame();
  
  console.log('RVM Visualization initialized successfully!');
}

// Start the application when the page loads
window.addEventListener('load', initialize);

// ===================================================================
// EMBEDDED WEBGL PLOTTING LIBRARY
// ===================================================================

(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
  typeof define === 'function' && define.amd ? define(['exports'], factory) :
  (global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.WebglPlotBundle = {}));
})(this, (function (exports) { 'use strict';

  class ColorRGBA {
      r;
      g;
      b;
      a;
      constructor(r, g, b, a) {
          this.r = r;
          this.g = g;
          this.b = b;
          this.a = a;
      }
  }

  /**
   * Baseline class
   */
  class WebglBase {
      //private static program: WebGLProgram;
      intensity;
      visible;
      /**
       * The number of data point pairs in the line
       */
      numPoints;
      /**
       * The data ponits for webgl array
       * @internal
       */
      xy;
      /**
       * The Color of the line
       */
      color;
      /**
       * The horizontal scale of the line
       * @default = 1
       */
      scaleX;
      /**
       * The vertical scale of the line
       * @default = 1
       */
      scaleY;
      /**
       * The horizontal offset of the line
       * @default = 0
       */
      offsetX;
      /**
       * the vertical offset of the line
       * @default = 0
       */
      offsetY;
      /**
       * if this is a close loop line or not
       * @default = false
       */
      loop;
      /**
       * total webgl number of points
       * @internal
       */
      webglNumPoints;
      /**
       * @private
       * @internal
       */
      _vbuffer;
      /**
       * @private
       * @internal
       */
      //public _prog: WebGLProgram;
      /**
       * @private
       * @internal
       */
      _coord;
      /**
       * @internal
       */
      constructor() {
          this.scaleX = 1;
          this.scaleY = 1;
          this.offsetX = 0;
          this.offsetY = 0;
          this.loop = false;
          this._vbuffer = 0;
          this._coord = 0;
          this.visible = true;
          this.intensity = 1;
          this.xy = new Float32Array([]);
          this.numPoints = 0;
          this.color = new ColorRGBA(0, 0, 0, 1);
          this.webglNumPoints = 0;
      }
  }

  /**
   * The standard Line class
   */
  class WebglLine extends WebglBase {
      currentIndex = 0;
      /**
       * Create a new line
       * @param c - the color of the line
       * @param numPoints - number of data pints
       * @example
       * ```typescript
       * x= [0,1]
       * y= [1,2]
       * line = new WebglLine( new ColorRGBA(0.1,0.1,0.1,1), 2);
       * ```
       */
      constructor(c, numPoints) {
          super();
          this.webglNumPoints = numPoints;
          this.numPoints = numPoints;
          this.color = c;
          this.xy = new Float32Array(2 * this.webglNumPoints);
      }
      /**
       * Set the X value at a specific index
       * @param index - the index of the data point
       * @param x - the horizontal value of the data point
       */
      setX(index, x) {
          this.xy[index * 2] = x;
      }
      /**
       * Set the Y value at a specific index
       * @param index : the index of the data point
       * @param y : the vertical value of the data point
       */
      setY(index, y) {
          this.xy[index * 2 + 1] = y;
      }
      /**
       * Get an X value at a specific index
       * @param index - the index of X
       */
      getX(index) {
          return this.xy[index * 2];
      }
      /**
       * Get an Y value at a specific index
       * @param index - the index of Y
       */
      getY(index) {
          return this.xy[index * 2 + 1];
      }
      /**
       * Make an equally spaced array of X points
       * @param start  - the start of the series
       * @param stepSize - step size between each data point
       *
       * @example
       * ```typescript
       * //x = [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]
       * const numX = 10;
       * line.lineSpaceX(-1, 2 / numX);
       * ```
       */
      lineSpaceX(start, stepSize) {
          for (let i = 0; i < this.numPoints; i++) {
              // set x to -num/2:1:+num/2
              this.setX(i, start + stepSize * i);
          }
      }
      /**
       * Automatically generate X between -1 and 1
       * equal to lineSpaceX(-1, 2/ number of points)
       */
      arrangeX() {
          this.lineSpaceX(-1, 2 / this.numPoints);
      }
      /**
       * Set a constant value for all Y values in the line
       * @param c - constant value
       */
      constY(c) {
          for (let i = 0; i < this.numPoints; i++) {
              // set x to -num/2:1:+num/2
              this.setY(i, c);
          }
      }
      /**
       * Add a new Y values to the end of current array and shift it, so that the total number of the pair remains the same
       * @param data - the Y array
       *
       * @example
       * ```typescript
       * yArray = new Float32Array([3, 4, 5]);
       * line.shiftAdd(yArray);
       * ```
       */
      shiftAdd(data) {
          const shiftSize = data.length;
          for (let i = 0; i < this.numPoints - shiftSize; i++) {
              this.setY(i, this.getY(i + shiftSize));
          }
          for (let i = 0; i < shiftSize; i++) {
              this.setY(i + this.numPoints - shiftSize, data[i]);
          }
      }
      /**
       * Add new Y values to the line and maintain the position of the last data point
       */
      addArrayY(yArray) {
          if (this.currentIndex + yArray.length <= this.numPoints) {
              for (let i = 0; i < yArray.length; i++) {
                  this.setY(this.currentIndex, yArray[i]);
                  this.currentIndex++;
              }
          }
      }
      /**
       * Replace the all Y values of the line
       */
      replaceArrayY(yArray) {
          if (yArray.length == this.numPoints) {
              for (let i = 0; i < this.numPoints; i++) {
                  this.setY(i, yArray[i]);
              }
          }
      }
  }

  /**
   * The step based line plot
   */
  class WebglStep extends WebglBase {
      /**
       * Create a new step line
       * @param c - the color of the line
       * @param numPoints - number of data pints
       * @example
       * ```typescript
       * x= [0,1]
       * y= [1,2]
       * line = new WebglStep( new ColorRGBA(0.1,0.1,0.1,1), 2);
       * ```
       */
      constructor(c, num) {
          super();
          this.webglNumPoints = num * 2;
          this.numPoints = num;
          this.color = c;
          this.xy = new Float32Array(2 * this.webglNumPoints);
      }
      /**
       * Set the Y value at a specific index
       * @param index - the index of the data point
       * @param y - the vertical value of the data point
       */
      setY(index, y) {
          this.xy[index * 4 + 1] = y;
          this.xy[index * 4 + 3] = y;
      }
      getX(index) {
          return this.xy[index * 4];
      }
      /**
       * Get an X value at a specific index
       * @param index - the index of X
       */
      getY(index) {
          return this.xy[index * 4 + 1];
      }
      /**
       * Make an equally spaced array of X points
       * @param start  - the start of the series
       * @param stepSize - step size between each data point
       *
       * @example
       * ```typescript
       * //x = [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]
       * const numX = 10;
       * line.lineSpaceX(-1, 2 / numX);
       * ```
       */
      lineSpaceX(start, stepsize) {
          for (let i = 0; i < this.numPoints; i++) {
              // set x to -num/2:1:+num/2
              this.xy[i * 4] = start + i * stepsize;
              this.xy[i * 4 + 2] = start + (i * stepsize + stepsize);
          }
      }
      /**
       * Set a constant value for all Y values in the line
       * @param c - constant value
       */
      constY(c) {
          for (let i = 0; i < this.numPoints; i++) {
              // set x to -num/2:1:+num/2
              this.setY(i, c);
          }
      }
      /**
       * Add a new Y values to the end of current array and shift it, so that the total number of the pair remains the same
       * @param data - the Y array
       *
       * @example
       * ```typescript
       * yArray = new Float32Array([3, 4, 5]);
       * line.shiftAdd(yArray);
       * ```
       */
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
      numPoints;
      xy;
      color;
      intenisty;
      visible;
      offsetTheta;
      constructor(c, numPoints) {
          super();
          this.webglNumPoints = numPoints;
          this.numPoints = numPoints;
          this.color = c;
          this.intenisty = 1;
          this.xy = new Float32Array(2 * this.webglNumPoints);
          this._vbuffer = 0;
          this._coord = 0;
          this.visible = true;
          this.offsetTheta = 0;
      }
      /**
       * @param index: index of the line
       * @param theta : angle in deg
       * @param r : radius
       */
      setRtheta(index, theta, r) {
          //const rA = Math.abs(r);
          //const thetaA = theta % 360;
          const x = r * Math.cos((2 * Math.PI * (theta + this.offsetTheta)) / 360);
          const y = r * Math.sin((2 * Math.PI * (theta + this.offsetTheta)) / 360);
          //const index = Math.round( ((theta % 360)/360) * this.numPoints );
          this.setX(index, x);
          this.setY(index, y);
      }
      getTheta(index) {
          //return Math.tan
          return 0;
      }
      getR(index) {
          //return Math.tan
          return Math.sqrt(Math.pow(this.getX(index), 2) + Math.pow(this.getY(index), 2));
      }
      setX(index, x) {
          this.xy[index * 2] = x;
      }
      setY(index, y) {
          this.xy[index * 2 + 1] = y;
      }
      getX(index) {
          return this.xy[index * 2];
      }
      getY(index) {
          return this.xy[index * 2 + 1];
      }
  }

  /**
   * The Square class
   */
  class WebglSquare extends WebglBase {
      /**
       * Create a new line
       * @param c - the color of the line
       * @example
       * ```typescript
       * line = new WebglSquare( new ColorRGBA(0.1,0.1,0.1,0.5) );
       * ```
       */
      constructor(c) {
          super();
          this.webglNumPoints = 4;
          this.numPoints = 4;
          this.color = c;
          this.xy = new Float32Array(2 * this.webglNumPoints);
      }
      /**
       * draw a square
       * @param x1 start x
       * @param y1 start y
       * @param x2 end x
       * @param y2 end y
       */
      setSquare(x1, y1, x2, y2) {
          this.xy = new Float32Array([x1, y1, x1, y2, x2, y1, x2, y2]);
      }
  }

  /**
   * modified functions from:
   * https://github.com/stackgl/gl-vec2
   * See License2.md for more info
   */
  const scaleAndAdd = (a, b, scale) => {
      const out = { x: 0, y: 0 };
      out.x = a.x + b.x * scale;
      out.y = a.y + b.y * scale;
      return out;
  };
  const normal = (dir) => {
      //get perpendicular
      const out = set(-dir.y, dir.x);
      return out;
  };
  const direction = (a, b) => {
      //get unit dir of two lines
      let out = subtract(a, b);
      out = normalize(out);
      return out;
  };
  /**
   * Adds two vec2's
   *
   * @param {vec2} out the receiving vector
   * @param {vec2} a the first operand
   * @param {vec2} b the second operand
   * @returns {vec2} out
   */
  const add = (a, b) => {
      const out = { x: 0, y: 0 };
      out.x = a.x + b.x;
      out.y = a.y + b.y;
      return out;
  };
  /**
   * Calculates the dot product of two vec2's
   *
   * @param {vec2} a the first operand
   * @param {vec2} b the second operand
   * @returns {Number} dot product of a and b
   */
  const dot = (a, b) => {
      return a.x * b.x + a.y * b.y;
  };
  /**
   * Normalize a vec2
   *
   * @param {vec2} out the receiving vector
   * @param {vec2} a vector to normalize
   * @returns {vec2} out
   */
  const normalize = (a) => {
      const out = { x: 0, y: 0 };
      let len = a.x * a.x + a.y * a.y;
      if (len > 0) {
          //TODO: evaluate use of glm_invsqrt here?
          len = 1 / Math.sqrt(len);
          out.x = a.x * len;
          out.y = a.y * len;
      }
      return out;
  };
  /**
   * Set the components of a vec2 to the given values
   *
   * @param {vec2} out the receiving vector
   * @param {Number} x X component
   * @param {Number} y Y component
   * @returns {vec2} out
   */
  const set = (x, y) => {
      const out = { x: 0, y: 0 };
      out.x = x;
      out.y = y;
      return out;
  };
  /**
   * Subtracts vector b from vector a
   *
   * @param {vec2} out the receiving vector
   * @param {vec2} a the first operand
   * @param {vec2} b the second operand
   * @returns {vec2} out
   */
  const subtract = (a, b) => {
      const out = { x: 0, y: 0 };
      out.x = a.x - b.x;
      out.y = a.y - b.y;
      return out;
  };

  /**
   * inspired and modified from:
   * https://github.com/mattdesl/polyline-normals
   * See License1.md for more info
   */
  const PolyLine = (lineXY) => {
      let curNormal;
      let lineA = { x: 0, y: 0 };
      let lineB = { x: 0, y: 0 };
      const out = [];
      const addNext = (normal, length) => {
          out.push({ vec2: normal, miterLength: length });
      };
      const getXY = (index) => {
          return { x: lineXY[index * 2], y: lineXY[index * 2 + 1] };
      };
      // add initial normals
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
          //stores tangent & miter
          const miter = computeMiter(lineA, lineB);
          const miterLen = computeMiterLen(lineA, miter, 1);
          addNext(miter, miterLen);
      }
      // add last normal
      // no miter, simple segment
      lineA = direction(getXY(numPoints - 1), getXY(numPoints - 2));
      curNormal = normal(lineA); //reset normal
      addNext(curNormal, 1);
      return out;
  };
  const computeMiter = (lineA, lineB) => {
      //get tangent line
      let tangent = add(lineA, lineB);
      tangent = normalize(tangent);
      //get miter as a unit vector
      const miter = set(-tangent.y, tangent.x);
      return miter;
  };
  const computeMiterLen = (lineA, miter, halfThick) => {
      const tmp = set(-lineA.y, lineA.x);
      //get the necessary length of our miter
      return halfThick / dot(miter, tmp);
  };

  /**
   * The standard Line class
   */
  class WebglThickLine extends WebglBase {
      currentIndex = 0;
      //protected triPoints: Float32Array;
      _linePoints;
      _thicknessRequested = 0;
      _actualThickness = 0;
      /**
       * Create a new line
       * @param c - the color of the line
       * @param numPoints - number of data pints
       * @example
       * ```typescript
       * x= [0,1]
       * y= [1,2]
       * line = new WebglLine( new ColorRGBA(0.1,0.1,0.1,1), 2);
       * ```
       */
      constructor(c, numPoints, thickness) {
          super();
          this.webglNumPoints = numPoints * 2;
          this.numPoints = numPoints;
          this.color = c;
          this._thicknessRequested = thickness;
          this._linePoints = new Float32Array(numPoints * 2);
          //this.triPoints = new Float32Array(this.numPoints * 4);
          this.xy = new Float32Array(2 * this.webglNumPoints);
      }
      convertToTriPoints() {
          //const thick = 0.01;
          const halfThick = this._actualThickness / 2;
          const normals = PolyLine(this._linePoints);
          //console.log(this.linePoints);
          //console.log(normals);
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
      /**
       * Set the X value at a specific index
       * @param index - the index of the data point
       * @param x - the horizontal value of the data point
       */
      setX(index, x) {
          this._linePoints[index * 2] = x;
      }
      /**
       * Set the Y value at a specific index
       * @param index : the index of the data point
       * @param y : the vertical value of the data point
       */
      setY(index, y) {
          this._linePoints[index * 2 + 1] = y;
      }
      /**
       * Make an equally spaced array of X points
       * @param start  - the start of the series
       * @param stepSize - step size between each data point
       *
       * @example
       * ```typescript
       * //x = [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]
       * const numX = 10;
       * line.lineSpaceX(-1, 2 / numX);
       * ```
       */
      lineSpaceX(start, stepSize) {
          for (let i = 0; i < this.numPoints; i++) {
              // set x to -num/2:1:+num/2
              this.setX(i, start + stepSize * i);
          }
      }
      setThickness(thickness) {
          this._thicknessRequested = thickness;
      }
      getThickness() {
          return this._thicknessRequested;
      }
      setActualThickness(thickness) {
          this._actualThickness = thickness;
      }
  }

  /**
   * Author Danial Chitnis 2019-20
   *
   * inspired by:
   * https://codepen.io/AzazelN28
   * https://www.tutorialspoint.com/webgl/webgl_modes_of_drawing.htm
   */
  /**
   * The main class for the webgl-plot library
   */
  class WebglPlot {
      /**
       * @private
       */
      webgl;
      /**
       * Global horizontal scale factor
       * @default = 1.0
       */
      gScaleX;
      /**
       * Global vertical scale factor
       * @default = 1.0
       */
      gScaleY;
      /**
       * Global X/Y scale ratio
       * @default = 1
       */
      gXYratio;
      /**
       * Global horizontal offset
       * @default = 0
       */
      gOffsetX;
      /**
       * Global vertical offset
       * @default = 0
       */
      gOffsetY;
      /**
       * Global log10 of x-axis
       * @default = false
       */
      gLog10X;
      /**
       * Global log10 of y-axis
       * @default = false
       */
      gLog10Y;
      /**
       * collection of data lines in the plot
       */
      _linesData;
      /**
       * collection of auxiliary lines (grids, markers, etc) in the plot
       */
      _linesAux;
      _thickLines;
      _surfaces;
      get linesData() {
          return this._linesData;
      }
      get linesAux() {
          return this._linesAux;
      }
      get thickLines() {
          return this._thickLines;
      }
      get surfaces() {
          return this._surfaces;
      }
      _progLine;
      /**
       * log debug output
       */
      debug = false;
      /**
       * Create a webgl-plot instance
       * @param canvas - the canvas in which the plot appears
       * @param debug - (Optional) log debug messages to console
       *
       * @example
       *
       * For HTMLCanvas
       * ```typescript
       * const canvas = document.getElementbyId("canvas");
       *
       * const devicePixelRatio = window.devicePixelRatio || 1;
       * canvas.width = canvas.clientWidth * devicePixelRatio;
       * canvas.height = canvas.clientHeight * devicePixelRatio;
       *
       * const webglp = new WebGLplot(canvas);
       * ...
       * ```
       * @example
       *
       * For OffScreenCanvas
       * ```typescript
       * const offscreen = htmlCanvas.transferControlToOffscreen();
       *
       * offscreen.width = htmlCanvas.clientWidth * window.devicePixelRatio;
       * offscreen.height = htmlCanvas.clientHeight * window.devicePixelRatio;
       *
       * const worker = new Worker("offScreenCanvas.js", { type: "module" });
       * worker.postMessage({ canvas: offscreen }, [offscreen]);
       * ```
       * Then in offScreenCanvas.js
       * ```typescript
       * onmessage = function (evt) {
       * const wglp = new WebGLplot(evt.data.canvas);
       * ...
       * }
       * ```
       */
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
          this._linesData = [];
          this._linesAux = [];
          this._thickLines = [];
          this._surfaces = [];
          //this.webgl = webgl;
          this.gScaleX = 1;
          this.gScaleY = 1;
          this.gXYratio = 1;
          this.gOffsetX = 0;
          this.gOffsetY = 0;
          this.gLog10X = false;
          this.gLog10Y = false;
          // Clear the color
          this.webgl.clear(this.webgl.COLOR_BUFFER_BIT);
          // Set the view port
          this.webgl.viewport(0, 0, canvas.width, canvas.height);
          this._progLine = this.webgl.createProgram();
          this.initThinLineProgram();
          //https://learnopengl.com/Advanced-OpenGL/Blending
          this.webgl.enable(this.webgl.BLEND);
          this.webgl.blendFunc(this.webgl.SRC_ALPHA, this.webgl.ONE_MINUS_SRC_ALPHA);
      }
      /**
       * updates and redraws the content of the plot
       */
      _drawLines(lines) {
          const webgl = this.webgl;
          lines.forEach((line) => {
              if (line.visible) {
                  webgl.useProgram(this._progLine);
                  const uscale = webgl.getUniformLocation(this._progLine, "uscale");
                  webgl.uniformMatrix2fv(uscale, false, new Float32Array([
                      line.scaleX * this.gScaleX * (this.gLog10X ? 1 / Math.log(10) : 1),
                      0,
                      0,
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
                      square.scaleX * this.gScaleX * (this.gLog10X ? 1 / Math.log(10) : 1),
                      0,
                      0,
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
              thickLine.scaleX * this.gScaleX * (this.gLog10X ? 1 / Math.log(10) : 1),
              0,
              0,
              thickLine.scaleY * this.gScaleY * this.gXYratio * (this.gLog10Y ? 1 / Math.log(10) : 1),
          ]));
          const uoffset = webgl.getUniformLocation(this._progLine, "uoffset");
          webgl.uniform2fv(uoffset, new Float32Array([thickLine.offsetX + this.gOffsetX, thickLine.offsetY + this.gOffsetY]));
          const isLog = webgl.getUniformLocation(this._progLine, "is_log");
          webgl.uniform2iv(isLog, new Int32Array([0, 0]));
          const uColor = webgl.getUniformLocation(this._progLine, "uColor");
          webgl.uniform4fv(uColor, [
              thickLine.color.r,
              thickLine.color.g,
              thickLine.color.b,
              thickLine.color.a,
          ]);
          webgl.drawArrays(webgl.TRIANGLE_STRIP, 0, thickLine.xy.length / 2);
      }
      _drawThickLines() {
          this._thickLines.forEach((thickLine) => {
              if (thickLine.visible) {
                  const calibFactor = Math.min(this.gScaleX, this.gScaleY);
                  //const calibFactor = 10;
                  //console.log(thickLine.getThickness());
                  thickLine.setActualThickness(thickLine.getThickness() / calibFactor);
                  thickLine.convertToTriPoints();
                  this._drawTriangles(thickLine);
              }
          });
      }
      /**
       * Draw and clear the canvas
       */
      update() {
          this.clear();
          this.draw();
      }
      /**
       * Draw without clearing the canvas
       */
      draw() {
          this._drawLines(this.linesData);
          this._drawLines(this.linesAux);
          this._drawThickLines();
          this._drawSurfaces(this.surfaces);
      }
      /**
       * Clear the canvas
       */
      clear() {
          //this.webgl.clearColor(0.1, 0.1, 0.1, 1.0);
          this.webgl.clear(this.webgl.COLOR_BUFFER_BIT);
      }
      /**
       * adds a line to the plot
       * @param line - this could be any of line, linestep, histogram, or polar
       *
       * @example
       * ```typescript
       * const line = new line(color, numPoints);
       * wglp.addLine(line);
       * ```
       */
      _addLine(line) {
          //line.initProgram(this.webgl);
          line._vbuffer = this.webgl.createBuffer();
          this.webgl.bindBuffer(this.webgl.ARRAY_BUFFER, line._vbuffer);
          this.webgl.bufferData(this.webgl.ARRAY_BUFFER, line.xy, this.webgl.STREAM_DRAW);
          //this.webgl.bindBuffer(this.webgl.ARRAY_BUFFER, line._vbuffer);
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
          // Create a vertex shader object
          const vertShader = this.webgl.createShader(this.webgl.VERTEX_SHADER);
          // Attach vertex shader source code
          this.webgl.shaderSource(vertShader, vertCode);
          // Compile the vertex shader
          this.webgl.compileShader(vertShader);
          // Fragment shader source code
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
      /**
       * remove the last data line
       */
      popDataLine() {
          this.linesData.pop();
      }
      /**
       * remove all the lines
       */
      removeAllLines() {
          this._linesData = [];
          this._linesAux = [];
          this._thickLines = [];
          this._surfaces = [];
      }
      /**
       * remove all data lines
       */
      removeDataLines() {
          this._linesData = [];
      }
      /**
       * remove all auxiliary lines
       */
      removeAuxLines() {
          this._linesAux = [];
      }
      /**
       * Change the WbGL viewport
       * @param a
       * @param b
       * @param c
       * @param d
       */
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