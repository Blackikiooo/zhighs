const std = @import("std");
const config = @import("config");

pub const HInt = switch (config.highs_int_width) {
    .w_32 => i32,
    .w_64 => i64,
};

pub const HUInt = switch (config.highs_int_width) {
    .w_32 => u32,
    .w_64 => u64,
};

test "public Highs integer types match configured width" {
    switch (config.highs_int_width) {
        .w_32 => {
            try std.testing.expect(HInt == i32);
            try std.testing.expect(HUInt == u32);
        },
        .w_64 => {
            try std.testing.expect(HInt == i64);
            try std.testing.expect(HUInt == u64);
        },
    }
}
