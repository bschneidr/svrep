as_data_frame_with_weights <- function(rep_design, full_weights,
                                       rep_wts_prefix) {

  # Extract weights from the desgin object
  replicate_weights <- stats::weights(weights, type = "analysis")
  sampling_weights <- stats::weights(weights, type = "sampling")


}
