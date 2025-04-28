

# SuperRiesz::super_riesz
# function (data, library, alternatives = list(), extra = list(), 
#           m = function(alpha, data) alpha(data()), folds = 5) 
# {
#     checkmate::assert_vector(library, min.len = 1)
#     checkmate::assert_function(m)
#     checkmate::assert_int(folds, lower = 1)
#     checkmate::assert_list(alternatives)
#     checkmate::assert_list(extra)
#     discrete <- TRUE
#     checkmate::assert_logical(discrete)
#     task <- TaskRiesz$new(id = "superriesz", backend = data, 
#                           alternatives = alternatives, extra = extra, m = m)
#     if (is.list(library)) {
#         learners <- lapply(library, function(learner) {
#             args <- learner[-1]
#             args$.key <- paste0("riesz.", learner[[1]])
#             do.call(lrn, args)
#         })
#     }
#     else {
#         learners <- lapply(library, function(learner) {
#             lrn(paste0("riesz.", learner))
#         })
#     }
#     if (is.null(folds)) 
#         folds = 5
#     if (folds > 1) {
#         cv <- riesz_superlearner_weights(learners, task, folds)
#         cv_risks <- cv$risks
#         cv_preds <- cv$preds
#     }
#     lapply(learners, function(learner) learner$train(task))
#     if (folds == 1) {
#         cv_risks <- unlist(lapply(learners, function(learner) learner$loss(task)))
#         names(cv_risks) <- unlist(lapply(learners, function(learner) learner$id))
#         cv_preds <- matrix(ncol = length(learners), nrow = task$nrow, 
#                            unlist(lapply(learners, function(learner) learner$predict(task)$response)))
#     }
#     if (length(library) == 1 || discrete) {
#         weights <- numeric(length(library))
#         weights[which.min(cv_risks)] <- 1
#     }
#     else {
#     }
#     sl <- list(learners = learners, weights = weights, m = m, 
#                risk = cv_risks)
#     class(sl) <- "super_riesz"
#     sl
# }

# nn_sequential_riesz_representer
# crumble::sequential_module
