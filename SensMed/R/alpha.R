
# adapted from crumble R package https://github.com/nt-williams/crumble/tree/main
hat.alpha <- function(train, valid, vars, architecture, .f, target, weights = NULL, control) {
    model <- nn_sequential_riesz_representer(
        train = train,
        x = vars,
        architecture = architecture,
        .f = .f,
        weights = weights,
        batch_size = control$batch_size,
        learning_rate = control$learning_rate,
        epochs = control$epochs,
        device = control$device
    )

    list(
        train = as.numeric(
            model(
                as_torch(
                    one_hot_encode(train[["data"]][, vars]),
                    device = control$device
                )
            )
        ),
        valid = as.numeric(
            model(
                as_torch(
                    one_hot_encode(valid[["data"]][, vars]),
                    device = control$device
                )
            )
        ),
        valid_target = as.numeric(
            model(
                as_torch(
                    one_hot_encode(valid[[ target ]][, vars]),
                    device = control$device
                )
            )
        )
    )

}

# from crumble R package https://github.com/nt-williams/crumble/tree/main
nn_sequential_riesz_representer <- function(train,
                                            x,
                                            architecture,
                                            .f,
                                            weights = NULL,
                                            batch_size,
                                            learning_rate,
                                            epochs,
                                            device) {
    dataset <- make_dataset(train, x, device = device)
    train_dl <- torch::dataloader(dataset, batch_size = batch_size)
    model <- architecture(ncol(dataset$data))
    model$to(device = device)

    weights <- weights %??% 1 # checkmate

    optimizer <- torch::optim_adam(
        params = c(model$parameters),
        lr = learning_rate,
        weight_decay = 0.01
    )

    scheduler <- torch::lr_one_cycle(
        optimizer,
        max_lr = learning_rate,
        total_steps = epochs
    )

    p <- progressr::progressor(steps = epochs)

    for (epoch in 1:epochs) {
        coro::loop(for (b in train_dl) {
            # Regression loss
            # loss <- (model(b$data)$pow(2) - (2 * weights * .f(model, b)))$mean(dtype = torch::torch_float())
            loss <- (model(b$data)$pow(2) - (2 * weights * .f(model, b)))$mean(dtype = torch::torch_float())

            optimizer$zero_grad()
            loss$backward()

            optimizer$step()
        })
        scheduler$step()
        p()
    }

    model$eval()
    model
}


make_dataset <- function(train, x, device = "cpu") {
    self <- NULL
    dataset <- torch::dataset(
        name = "tmp_dataset",
        initialize = function(train, x, device) {
            for (df in names(train)) {
                if (ncol(train[[df]]) > 0) {
                    df_x <- train[[df]][, x, drop = FALSE]
                    # self[[df]] <- torch::torch_tensor(as.matrix(df_x), dtype = torch::torch_float(), device = device)
                    self[[df]] <- one_hot_encode(df_x) |>
                        as_torch(device = device)
                }
            }
        },
        .getitem = function(i) {
            fields <- grep("data", names(self), value = TRUE)
            setNames(lapply(fields, function(x) self[[x]][i, ]), fields)
        },
        .length = function() {
            self$data$size()[1]
        }
    )
    dataset(train, x, device)
}


sequential_module <- function(layers = 1, hidden = 20, dropout = 0.1) {
    function(d_in) {
        d_out <- 1

        middle_layers <- lapply(1:layers, \(x) torch::nn_sequential(torch::nn_linear(hidden, hidden), torch::nn_elu()))

        torch::nn_sequential(
            torch::nn_linear(d_in, hidden),
            torch::nn_elu(),
            do.call(torch::nn_sequential, middle_layers),
            torch::nn_linear(hidden, d_out),
            torch::nn_dropout(dropout),
            torch::nn_softplus()
        )
    }
}
