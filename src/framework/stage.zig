//! Solver lifecycle stages shared by API validation and component scheduling.

pub const Stage = enum {
    empty,
    problem,
    transformed,
    presolving,
    presolved,
    solving,
    solved,
    deinitialized,
};
