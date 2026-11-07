# The blrm prediction and recursion code, checked draw by draw against the readable
# loops in helper-slow-loop.R. Sampling from Stan is too slow for the default run,
# so these fits are made by hand: the draws are real numbers, the transition
# probabilities depend on time, treatment and both previous states, and only the
# design matrix is mocked. test-slow-blrm.R covers the design matrix with real fits.

markov_fake_blrm <- function(partial_odds = FALSE) {
  draws <- rbind(
    c(1.2, 0.1, -1.3, 0.40, 0.08, -0.60, -0.25),
    c(0.9, -0.2, -1.1, 0.10, 0.03, -0.45, -0.35),
    c(1.5, 0.3, -1.6, 0.70, 0.12, -0.75, -0.15),
    c(1.1, 0.0, -1.2, 0.25, 0.05, -0.55, -0.30)
  )
  colnames(draws) <- c("y>=2", "y>=3", "y>=4", "tx", "time", "yprev", "ypprev")
  if (partial_odds) {
    draws <- cbind(draws, tau_time = c(0.05, -0.10, 0.15, 0.02))
  }

  structure(
    list(
      draws = draws,
      non.slopes = 3L,
      pppo = as.integer(partial_odds),
      tauInfo = data.frame(name = if (partial_odds) "tau_time" else character()),
      cppo = if (partial_odds) function(y) as.numeric(y) - 1,
      ylevels = 1:4,
      yname = "y",
      clusterInfo = list(name = "id", cluster = c("a", "b"))
    ),
    class = c("blrm", "orm"),
    markov_engine = "blrm"
  )
}

markov_fake_design <- function(model, newdata, second = FALSE) {
  time <- as.numeric(as.character(newdata$time))
  if (second) {
    return(matrix(time, ncol = 1, dimnames = list(NULL, "tau_time")))
  }
  cbind(
    tx = newdata$tx,
    time = time,
    yprev = as.numeric(as.character(newdata$yprev)),
    ypprev = as.numeric(as.character(newdata$ypprev))
  )
}

# One draw's transition probabilities for one row, written from the model:
# P(Y >= k) = plogis(alpha_k + x * beta + time * tau * cppo(k)).
markov_fake_transition <- function(model, draw) {
  coefficient <- model$draws[draw, ]
  function(row) {
    x <- markov_fake_design(model, row)
    linear <- sum(x[1, ] * coefficient[colnames(x)])
    tau <- if (model$pppo > 0L) {
      row$time * coefficient[["tau_time"]] * model$cppo(2:4)
    } else {
      0
    }
    at_least <- c(1, stats::plogis(coefficient[1:3] + linear + tau), 0)
    -diff(at_least)
  }
}

markov_fake_profiles <- function() {
  data.frame(
    id = c("a", "b", "c", "d"),
    time = 1,
    tx = c(0, 1, 0, 1),
    marker = c(1L, 1L, 1L, 0L),
    yprev = factor(c(1, 2, 3, 2), levels = 1:4),
    ypprev = factor(c(3, 1, 2, 2), levels = 1:4)
  )
}

markov_fake_reference <- function(model, profiles, times, absorb, second_order) {
  n_draws <- nrow(model$draws)
  lapply(seq_len(n_draws), function(draw) {
    predict_transition <- markov_fake_transition(model, draw)
    lapply(seq_len(nrow(profiles)), function(i) {
      patient <- profiles[i, , drop = FALSE]
      if (second_order) {
        slow_second_order_sops(predict_transition, patient, times, 1:4, absorb)
      } else {
        slow_first_order_sops(predict_transition, patient, times, 1:4, absorb)
      }
    })
  })
}

markov_fake_cases <- list(
  list(label = "first order", partial_odds = FALSE, second_order = FALSE),
  list(label = "first order, partial odds", partial_odds = TRUE, second_order = FALSE),
  list(label = "second order", partial_odds = FALSE, second_order = TRUE),
  list(label = "second order, partial odds", partial_odds = TRUE, second_order = TRUE)
)

test_that("blrm sops() draws match the readable loop draw by draw", {
  local_mocked_bindings(blrm_design_matrix = markov_fake_design)
  profiles <- markov_fake_profiles()
  times <- 1:4

  for (case in markov_fake_cases) {
    model <- markov_fake_blrm(case$partial_odds)
    for (absorb in list(4L, c(3L, 4L))) {
      keep <- !as.character(profiles$yprev) %in% as.character(absorb)
      kept <- profiles[keep, , drop = FALSE]
      rownames(kept) <- NULL
      args <- list(
        model,
        newdata = kept,
        times = times,
        absorb = absorb,
        n_draws = NULL,
        return_draws = TRUE,
        posterior_summary = "mean"
      )
      if (case$second_order) {
        args$p2_var <- "ypprev"
      }
      result <- do.call(sops, args)
      draws <- attr(result, "draws")
      reference <- markov_fake_reference(
        model,
        kept,
        times,
        absorb,
        case$second_order
      )

      expect_identical(
        nrow(draws),
        nrow(model$draws) * nrow(kept) * length(times) * 4L
      )
      gaps <- vapply(
        seq_len(nrow(draws)),
        function(row) {
          loop <- reference[[draws$draw_id[row]]][[draws$rowid[row]]]
          abs(
            draws$estimate[row] -
              loop[draws$time[row], as.integer(as.character(draws$state[row]))]
          )
        },
        numeric(1)
      )
      expect_lt(
        max(gaps),
        1e-9,
        label = paste(case$label, "absorbing", toString(absorb))
      )

      # The reported estimate is the mean over draws.
      means <- vapply(
        seq_len(nrow(result)),
        function(row) {
          mean(vapply(
            reference,
            function(one_draw) {
              one_draw[[result$rowid[row]]][
                result$time[row],
                as.integer(as.character(result$state[row]))
              ]
            },
            numeric(1)
          ))
        },
        numeric(1)
      )
      expect_equal(result$estimate, means, tolerance = 1e-10)
    }
  }
})

test_that("blrm avg_sops() draws match averaging the loop within unequal groups", {
  local_mocked_bindings(blrm_design_matrix = markov_fake_design)
  profiles <- markov_fake_profiles()
  times <- 1:3

  for (case in markov_fake_cases) {
    model <- markov_fake_blrm(case$partial_odds)
    args <- list(
      model,
      newdata = profiles,
      by = "marker",
      times = times,
      absorb = 4L,
      n_draws = NULL,
      return_draws = TRUE,
      posterior_summary = "mean"
    )
    if (case$second_order) {
      args$p2_var <- "ypprev"
    }
    result <- do.call(avg_sops, args)
    draws <- attr(result, "draws")
    reference <- markov_fake_reference(
      model,
      profiles,
      times,
      4L,
      case$second_order
    )

    gaps <- vapply(
      seq_len(nrow(draws)),
      function(row) {
        members <- which(profiles$marker == draws$marker[row])
        values <- vapply(
          reference[[draws$draw_id[row]]][members],
          function(loop) {
            loop[draws$time[row], as.integer(as.character(draws$state[row]))]
          },
          numeric(1)
        )
        abs(draws$estimate[row] - mean(values))
      },
      numeric(1)
    )
    expect_lt(max(gaps), 1e-9, label = case$label)
  }
})
