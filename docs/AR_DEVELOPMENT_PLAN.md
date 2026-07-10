# Oviz AR Development Plan

## Current status

AR is experimental and disabled by default. Standard mobile and desktop viewers do not include the AR button or exporter runtime. Development builds can opt in with:

```python
threejs_initial_state={"ar_enabled": True}
```

The static USDZ exporter remains available for experiments, but it is not the target for high-fidelity Edenhofer dust because USDZ cannot preserve Oviz's live volume ray-marching shader.

## Track 1: Camera-backed web mode

This is the fastest shareable path. It keeps Oviz's existing Three.js `Data3DTexture` and volume shader, then places the transparent WebGL canvas over the iPhone camera feed.

Planned work:

1. Add an experimental camera mode using `getUserMedia`.
2. Keep the current Oviz scene, lasso state, cluster points, and full-resolution volume renderer.
3. Use device orientation for camera rotation and touch gestures for placement, scale, and recentering.
4. Add clear camera-permission, unsupported-browser, and exit states.
5. Test rendering quality, memory use, thermal behavior, and interaction on a physical iPhone.

This mode can closely match Oviz's volume appearance and remains a normal shareable web page. It will not provide reliable ARKit plane anchors, room occlusion, or full positional tracking.

## Track 2: Native ARKit, then App Clip

This is the path to true anchored AR with the best achievable volumetric fidelity. Start with a normal development app; only package it as an App Clip after the renderer works well on-device.

Planned work:

1. Create a small Swift/ARKit test app and place a basic Oviz scene on a detected plane.
2. Port the Oviz volume ray marcher from GLSL to Metal rather than converting the volume to USDZ geometry.
3. Load and decompress the Edenhofer volume as an 8-bit 3D texture, initially at the viewer's current resolution.
4. Add cluster points, selection input, scale, placement, recentering, and scene reset.
5. Measure memory, frame rate, startup time, and thermal limits on a physical iPhone; add device-specific downsampling only if measurements require it.
6. Add the sky-dome mode after the 3D volume path is stable.
7. Create the parent iOS app and App Clip, configure the associated domain and invocation URL, then test distribution and review constraints.

Codex can create and maintain the Xcode project, Swift and Metal code, data conversion, builds, tests, and device-log analysis. A person must handle Apple account credentials and agreements, signing access, device trust and Developer Mode, camera permissions, physical movement during testing, and App Store Connect submission decisions.

## Decision gates

- Choose the camera-backed web mode when fidelity and link sharing matter more than physical anchoring.
- Choose native ARKit when stable plane placement, positional tracking, and room-aware behavior are required.
- Do not begin App Clip packaging until the normal ARKit app renders the Edenhofer volume acceptably on the target iPhone.
- Keep AR disabled in production Oviz exports until one path passes physical-device acceptance testing.
