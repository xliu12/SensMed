med_data <- S7::new_class("med_data",
                              properties = list(
                                  data = S7::new_property(S7::class_data.frame),
                                  vars = S7::new_property(S7::new_class("med_vars")),
                                  weights = S7::new_property(S7::class_numeric),
                                  a1 = S7::new_property(S7::class_numeric),
                                  a0 = S7::new_property(S7::class_numeric),
                                  # d0 = S7::new_property(S7::class_function, default = NULL),
                                  # d1 = S7::new_property(S7::class_function, default = NULL),
                                  data_0 = S7::new_property(S7::class_data.frame),
                                  data_1 = S7::new_property(S7::class_data.frame),
                                  data_indepM = S7::new_property(S7::class_data.frame),
                                  data_0_indepM = S7::new_property(S7::class_data.frame),
                                  data_1_indepM = S7::new_property(S7::class_data.frame)
                                  # data_0zp = S7::new_property(S7::class_data.frame),
                                  # data_1zp = S7::new_property(S7::class_data.frame)
                              ),
                              constructor = function(data, vars, weights, a1, a0) {
                                  
                                  S7::new_object(
                                      S7::S7_object(),
                                      data = data,
                                      vars = vars,
                                      weights = normalize(weights),
                                      id = 1:nrow(data),
                                      a1 = a1,
                                      a0 = a0,
                                      data_1 = shift_data(data, vars@A, a1, center = vars@S),
                                      data_0 = shift_data(data, vars@A, a0, center = vars@S),
                                      data_indepM = data,
                                      data_0_indepM = data,
                                      data_1_indepM = data
                                      # data_0zp = data.frame(),
                                      # data_1zp = data.frame()
                                  )
                              }
                              
)


med_vars <- S7::new_class("med_vars",
                              properties = list(
                                  A = S7::class_character,
                                  Y = S7::class_character,
                                  M = S7::class_character,
                                  C = S7::class_character,
                                  cluster_id = S7::new_property(class = S7::class_character, default = NA_character_),
                                  Wj = S7::new_property(class = S7::class_character, default = NA_character_),
                                  S = S7::new_property(class = S7::class_character, default = NA_character_),
                                  id = S7::new_property(class = S7::class_character, default = NA_character_)
                              )
)


training <- S7::new_generic("training", "x")
validation <- S7::new_generic("validation", "x")

S7::method(training, med_data) <- function(x, fold_obj, fold) {
    list(
        data = x@data[fold_obj[[fold]]$training_set, , drop = FALSE],
        data_0 = x@data_0[fold_obj[[fold]]$training_set, , drop = FALSE],
        data_1 = x@data_1[fold_obj[[fold]]$training_set, , drop = FALSE] 
        ,
        data_indepM = x@data_indepM[fold_obj[[fold]]$training_set, , drop = FALSE],
        data_0_indepM = x@data_0_indepM[fold_obj[[fold]]$training_set, , drop = FALSE],
        data_1_indepM = x@data_1_indepM[fold_obj[[fold]]$training_set, , drop = FALSE]
        # data_0zp = x@data_0zp[fold_obj[[fold]]$training_set, , drop = FALSE],
        # data_1zp = x@data_1zp[fold_obj[[fold]]$training_set, , drop = FALSE]
    )
}

S7::method(validation, med_data) <- function(x, fold_obj, fold) {
    list(
        data = x@data[fold_obj[[fold]]$validation_set, , drop = FALSE],
        data_0 = x@data_0[fold_obj[[fold]]$validation_set, , drop = FALSE],
        data_1 = x@data_1[fold_obj[[fold]]$validation_set, , drop = FALSE] 
        ,
        data_indepM = x@data_indepM[fold_obj[[fold]]$validation_set, , drop = FALSE],
        data_0_indepM = x@data_0_indepM[fold_obj[[fold]]$validation_set, , drop = FALSE],
        data_1_indepM = x@data_1_indepM[fold_obj[[fold]]$validation_set, , drop = FALSE]
        # data_0zp = x@data_0zp[fold_obj[[fold]]$validation_set, , drop = FALSE],
        # data_1zp = x@data_1zp[fold_obj[[fold]]$validation_set, , drop = FALSE]
    )
}
