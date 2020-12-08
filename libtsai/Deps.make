tf_calib_data.ho: tf_calib_data.h
tf_calib_data.ho: /home/jstolfi/include/r2.ho
tf_calib_data.ho: /home/jstolfi/include/r3.ho

tf_lmdif.ho: tf_lmdif.h

tf_errors.ho: tf_errors.h
tf_errors.ho: /home/jstolfi/include/r3.ho
tf_errors.ho: /home/jstolfi/include/r2.ho
tf_errors.ho: /home/jstolfi/include/r4x4.ho
tf_errors.ho: tf_camera.ho
tf_errors.ho: tf_calib.ho

tf_camera_specs.ho: tf_camera_specs.h
tf_camera_specs.ho: /home/jstolfi/include/r2.ho
tf_camera_specs.ho: /home/jstolfi/include/r3.ho
tf_camera_specs.ho: /home/jstolfi/include/r3x3.ho
tf_camera_specs.ho: /home/jstolfi/include/r4.ho
tf_camera_specs.ho: /home/jstolfi/include/r4x4.ho
tf_camera_specs.ho: /home/jstolfi/include/interval.ho
tf_camera_specs.ho: tf_camera.ho
tf_camera_specs.ho: tf_camera_specs.ho

tf_camera_plot.ho: tf_camera_plot.h
tf_camera_plot.ho: /home/jstolfi/include/pswr.ho
tf_camera_plot.ho: /home/jstolfi/include/frgb.ho
tf_camera_plot.ho: tf_camera.ho
tf_camera_plot.ho: tf_errors.ho
tf_camera_plot.ho: tf_targets.ho
tf_camera_plot.ho: tf_calib.ho

tf_calib_guess2.ho: tf_calib_guess2.h
tf_calib_guess2.ho: /home/jstolfi/include/r3.ho
tf_calib_guess2.ho: /home/jstolfi/include/r2.ho
tf_calib_guess2.ho: /home/jstolfi/include/r4x4.ho
tf_calib_guess2.ho: tf_camera.ho
tf_calib_guess2.ho: tf_calib.ho
tf_calib_guess2.ho: tf_calib_guess1.ho

tf_math.ho: tf_math.h

tf_calib.ho: tf_calib.h
tf_calib.ho: /home/jstolfi/include/r3.ho
tf_calib.ho: /home/jstolfi/include/r2.ho
tf_calib.ho: /home/jstolfi/include/r4x4.ho
tf_calib.ho: tf_camera.ho
tf_calib.ho: tf_camera_specs.ho
tf_calib.ho: tf_calib_data.ho

tf_camera_image.ho: tf_camera_image.h
tf_camera_image.ho: /home/jstolfi/include/float_image.ho
tf_camera_image.ho: tf_camera.ho

tf_camera.ho: tf_camera.h
tf_camera.ho: /home/jstolfi/include/r2.ho
tf_camera.ho: /home/jstolfi/include/r3.ho
tf_camera.ho: /home/jstolfi/include/r3x3.ho
tf_camera.ho: /home/jstolfi/include/r4.ho
tf_camera.ho: /home/jstolfi/include/r4x4.ho

tf_calib_guess1.ho: tf_calib_guess1.h
tf_calib_guess1.ho: /home/jstolfi/include/r3.ho
tf_calib_guess1.ho: /home/jstolfi/include/r2.ho
tf_calib_guess1.ho: /home/jstolfi/include/r4x4.ho
tf_calib_guess1.ho: /home/jstolfi/include/rmxn.ho
tf_calib_guess1.ho: tf_camera.ho
tf_calib_guess1.ho: tf_camera_specs.ho
tf_calib_guess1.ho: tf_matrix.ho
tf_calib_guess1.ho: tf_calib.ho

tf_calib_quality.ho: tf_calib_quality.h
tf_calib_quality.ho: tf_camera.ho
tf_calib_quality.ho: tf_calib.ho

tf_matrix.ho: tf_matrix.h
tf_matrix.ho: /home/jstolfi/include/r3.ho
tf_matrix.ho: /home/jstolfi/include/r2.ho
tf_matrix.ho: /home/jstolfi/include/r4x4.ho

tf_calib_refine2.ho: tf_calib_refine2.h
tf_calib_refine2.ho: /home/jstolfi/include/r3.ho
tf_calib_refine2.ho: /home/jstolfi/include/r2.ho
tf_calib_refine2.ho: /home/jstolfi/include/r4x4.ho
tf_calib_refine2.ho: tf_camera.ho
tf_calib_refine2.ho: tf_calib.ho

tf_kalman.ho: tf_kalman.h

tf_targets.ho: tf_targets.h
tf_targets.ho: /home/jstolfi/include/r2.ho
tf_targets.ho: /home/jstolfi/include/r3.ho
tf_targets.ho: tf_camera.ho
tf_targets.ho: tf_calib_data.ho

tf_targets.o: tf_targets.c
tf_targets.o: /home/jstolfi/include/r3.ho
tf_targets.o: /home/jstolfi/include/r2.ho
tf_targets.o: /home/jstolfi/include/affirm.ho
tf_targets.o: /home/jstolfi/include/jsfile.ho
tf_targets.o: tf_camera.ho
tf_targets.o: tf_calib.ho
tf_targets.o: tf_targets.ho
tf_targets.o: tf_math.ho

tf_lmdif.o: tf_lmdif.c
tf_lmdif.o: tf_lmdif.ho

tf_calib_data.o: tf_calib_data.c
tf_calib_data.o: /home/jstolfi/include/jsfile.ho
tf_calib_data.o: /home/jstolfi/include/affirm.ho
tf_calib_data.o: /home/jstolfi/include/r2.ho
tf_calib_data.o: /home/jstolfi/include/r3.ho
tf_calib_data.o: /home/jstolfi/include/fget.ho
tf_calib_data.o: tf_calib_data.ho

tf_camera_image.o: tf_camera_image.c
tf_camera_image.o: /home/jstolfi/include/float_image.ho
tf_camera_image.o: /home/jstolfi/include/float_image_interpolate.ho
tf_camera_image.o: /home/jstolfi/include/affirm.ho
tf_camera_image.o: tf_camera.ho
tf_camera_image.o: tf_camera_image.ho

tf_calib_guess1.o: tf_calib_guess1.c
tf_calib_guess1.o: /home/jstolfi/include/jsfile.ho
tf_calib_guess1.o: tf_camera.ho
tf_calib_guess1.o: tf_matrix.ho
tf_calib_guess1.o: tf_errors.ho
tf_calib_guess1.o: tf_math.ho
tf_calib_guess1.o: tf_calib.ho
tf_calib_guess1.o: tf_calib_guess1.ho
tf_calib_guess1.o: tf_calib_refine2.ho

tf_errors.o: tf_errors.c
tf_errors.o: tf_calib.ho
tf_errors.o: tf_camera.ho
tf_errors.o: tf_errors.ho
tf_errors.o: tf_math.ho

tf_camera_specs.o: tf_camera_specs.c
tf_camera_specs.o: /home/jstolfi/include/jsfile.ho
tf_camera_specs.o: /home/jstolfi/include/interval.ho
tf_camera_specs.o: /home/jstolfi/include/affirm.ho
tf_camera_specs.o: /home/jstolfi/include/r3.ho
tf_camera_specs.o: /home/jstolfi/include/fget.ho
tf_camera_specs.o: /home/jstolfi/include/r2.ho
tf_camera_specs.o: /home/jstolfi/include/r4x4.ho
tf_camera_specs.o: tf_camera.ho
tf_camera_specs.o: tf_camera_specs.ho

tf_calib_guess2.o: tf_calib_guess2.c
tf_calib_guess2.o: /home/jstolfi/include/jsfile.ho
tf_calib_guess2.o: /home/jstolfi/include/affirm.ho
tf_calib_guess2.o: tf_camera.ho
tf_calib_guess2.o: tf_matrix.ho
tf_calib_guess2.o: tf_errors.ho
tf_calib_guess2.o: tf_math.ho
tf_calib_guess2.o: tf_calib.ho
tf_calib_guess2.o: tf_calib_guess2.ho

tf_calib.o: tf_calib.c
tf_calib.o: /home/jstolfi/include/jsfile.ho
tf_calib.o: /home/jstolfi/include/affirm.ho
tf_calib.o: tf_lmdif.ho
tf_calib.o: tf_camera.ho
tf_calib.o: tf_calib.ho
tf_calib.o: tf_matrix.ho
tf_calib.o: tf_math.ho
tf_calib.o: tf_errors.ho

tf_camera.o: tf_camera.c
tf_camera.o: /home/jstolfi/include/jsfile.ho
tf_camera.o: /home/jstolfi/include/affirm.ho
tf_camera.o: /home/jstolfi/include/r3.ho
tf_camera.o: /home/jstolfi/include/r2.ho
tf_camera.o: /home/jstolfi/include/r4x4.ho
tf_camera.o: tf_camera.ho

tf_camera_plot.o: tf_camera_plot.c
tf_camera_plot.o: /home/jstolfi/include/affirm.ho
tf_camera_plot.o: /home/jstolfi/include/pswr.ho
tf_camera_plot.o: /home/jstolfi/include/jsfile.ho
tf_camera_plot.o: tf_camera.ho
tf_camera_plot.o: tf_errors.ho
tf_camera_plot.o: tf_targets.ho
tf_camera_plot.o: tf_calib.ho
tf_camera_plot.o: tf_camera_plot.ho

tf_matrix.o: tf_matrix.c
tf_matrix.o: /home/jstolfi/include/rmxn.ho
tf_matrix.o: tf_matrix.ho

tf_calib_refine2.o: tf_calib_refine2.c
tf_calib_refine2.o: /home/jstolfi/include/jsfile.ho
tf_calib_refine2.o: /home/jstolfi/include/affirm.ho
tf_calib_refine2.o: tf_camera.ho
tf_calib_refine2.o: tf_calib.ho
tf_calib_refine2.o: tf_matrix.ho
tf_calib_refine2.o: tf_math.ho
tf_calib_refine2.o: tf_calib.ho
tf_calib_refine2.o: tf_calib_refine2.ho

tf_kalman.o: tf_kalman.c
tf_kalman.o: /home/jstolfi/include/affirm.ho
tf_kalman.o: tf_kalman.ho

