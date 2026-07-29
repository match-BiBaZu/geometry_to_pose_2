"""
quat_probe.py

Single-trial quaternion probe for calibrating chute_drop_batch.py's pose
matching against ChutePoseAnalysis.m's refQuat table.

Drops ONE part into the chute at a given (or random) initial orientation,
waits for it to settle, opens the GUI so you can SEE which pose it landed
in, and prints:
  - the final orientation in WORLD frame
  - the final orientation in the CHUTE-LOCAL frame (same transform
    run_condition_trials uses before matching against refQuat)

You look at the GUI, find the matching row in the MASTER_summary /
GEOMETRIC_summary table by eye, and report back which Pose number it was.
That refQuat can then be compared against the printed local-frame
quaternion to check for a systematic frame/convention offset.

Reuses chute_drop_batch.py's own functions directly (mesh loading, VHACD,
COM computation, world setup, settle loop, frame transforms) so this is
guaranteed to use the exact same pipeline as the real batch run -- no
reimplementation, no chance of drifting out of sync with it.

--------------------------------------------------------------------------
USAGE (edit the CONFIG block below, then run)
--------------------------------------------------------------------------
    python quat_probe.py

Or from a notebook:
    import quat_probe
    quat_probe.run()
"""

import numpy as np
import pybullet as p

import chute_drop_batch as cdb


# ===========================================================================
# CONFIG -- edit these
# ===========================================================================

PART_MESH_PATH   = r"C:\Users\benny\OneDrive\Documents\LLM\PyBullet\SIM\Rk2i.stl"
CHUTE_PATH       = r"C:\Users\benny\OneDrive\Documents\LLM\PyBullet\SIM\chute.obj"
PART_SCALE       = 0.015
CHUTE_SCALE      = 1.0
CHUTE_CONCAVE    = True

ROLL_DEG  = 15.0
PITCH_DEG = 15.0

DROP_HEIGHT = 0.3          # meters, above chute origin along chute local +Z
DROP_XY     = [0.0, 0.0]   # meters, chute-local X/Y offset
PART_MASS   = 0.01

# Initial orientation dropped from -- None = uniformly random each run,
# or pass an explicit (x,y,z,w) tuple to always start from the same
# orientation (useful for reproducing a specific outcome).
INIT_ORN_XYZW = None

SEED       = None      # None = different every run
MAX_STEPS  = 20000
SLEEP      = 0.01       # slow down GUI playback so you can actually watch it settle

VHACD_CACHE_DIR = "quat_probe_vhacd_cache"


def run():
    rng = np.random.default_rng(SEED)

    p.connect(p.GUI)

    # --- Build the tilted chute world (reuses batch harness logic exactly) ---
    chute_id, chute_quat_xyzw = cdb.reset_condition_world(
        CHUTE_PATH, [CHUTE_SCALE] * 3, CHUTE_CONCAVE,
        ROLL_DEG, PITCH_DEG, catch_plane=True,
        z_lift=2.0, gui=True,
    )
    p.changeDynamics(chute_id, -1, lateralFriction=0.08, restitution=0.02,
                      rollingFriction=0.005, spinningFriction=0.005,
                      linearDamping=0.04, angularDamping=0.1)

    # --- Part mesh: convex-decompose + true center of mass, same as batch ---
    part_collision_path = cdb.get_convex_decomposition(PART_MESH_PATH, VHACD_CACHE_DIR)
    part_com = cdb.compute_part_center_of_mass(PART_MESH_PATH, PART_SCALE)

    # --- Drop position: chute-local offset, rotated + lifted like the batch script ---
    drop_local_offset = [DROP_XY[0], DROP_XY[1], DROP_HEIGHT]
    drop_pos_offset = cdb.rotate_vec_by_quat_xyzw(drop_local_offset, chute_quat_xyzw)
    drop_pos = [drop_pos_offset[0], drop_pos_offset[1], drop_pos_offset[2] + 2.0]

    orn0 = INIT_ORN_XYZW if INIT_ORN_XYZW is not None else cdb.random_quaternion_xyzw(rng)

    body_id = cdb.spawn_part(part_collision_path, [PART_SCALE] * 3, PART_MASS,
                              drop_pos, orn0, part_com=part_com)
    p.changeDynamics(body_id, -1, lateralFriction=0.08, restitution=0.02,
                      rollingFriction=0.005, spinningFriction=0.005,
                      linearDamping=0.05, angularDamping=0.1)

    print(f"\nDropping at roll={ROLL_DEG} pitch={PITCH_DEG}, init orn (xyzw) = {orn0}")
    print("Watching until it settles ...\n")

    steps, settled = cdb.run_until_settled(
        body_id, max_steps=MAX_STEPS, gui=True, sleep=SLEEP,
        ambient_jitter=True, weight=PART_MASS * cdb.GRAVITY_MAG,
        jitter_interval=200, jitter_force_frac=0.02, jitter_torque_frac=0.02,
        rng=rng)

    still_settled = True
    if settled:
        still_settled = cdb.verify_settled_not_knife_edge(
            body_id, mass=PART_MASS, rng=rng,
            nudge_force_frac=0.3, nudge_torque_frac=0.3,
            nudge_steps=10, recheck_steps=200, gui=True, sleep=SLEEP)

    final_pos, final_orn_world = p.getBasePositionAndOrientation(body_id)
    final_orn_local = cdb.rotate_to_chute_local(final_orn_world, chute_quat_xyzw)

    x, y, z, w = final_orn_local
    local_wxyz = (w, x, y, z)

    print("=" * 70)
    print(f"settled           : {settled}")
    print(f"knife-edge check  : {'passed' if still_settled else 'FAILED (still moving after nudge)'}")
    print(f"steps to settle   : {steps}")
    print(f"final world pos   : {final_pos}")
    print(f"final world orn (xyzw)  : {final_orn_world}")
    print(f"final LOCAL orn (xyzw)  : {final_orn_local}")
    print(f"final LOCAL orn (wxyz)  : ({local_wxyz[0]:.4f}, {local_wxyz[1]:.4f}, "
          f"{local_wxyz[2]:.4f}, {local_wxyz[3]:.4f})")
    print("=" * 70)
    print("\nCompare the (w,x,y,z) line above against the refQuat column in")
    print(f"Rk2i_{cdb.matlab_round(ROLL_DEG)}R_{cdb.matlab_round(PITCH_DEG)}P_MASTER_summary.txt / "
          f"_GEOMETRIC_summary.txt")
    print("Look at the GUI to see which pose it visually landed in, then tell")
    print("me the Pose number -- we'll check the geodesic distance and delta")
    print("rotation between the two.\n")
    print("(GUI window will stay open -- close it or Ctrl+C when done looking.)")

    while p.isConnected():
        p.stepSimulation()


if __name__ == "__main__":
    run()
