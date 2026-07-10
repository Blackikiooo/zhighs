//! Component categories understood by the registry and scheduler.

pub const Kind = enum {
    benders,
    benders_cut,
    constraint_handler,
    presolver,
    separator,
    propagator,
    branching_rule,
    primal_heuristic,
    variable_pricer,
    conflict_handler,
    iis_finder,
    node_selector,
    cut_selector,
    event_handler,
    reader,
    relaxation_handler,
};
