const std = @import("std");
const config = @import("config");

pub const HInt = switch (config.highs_int_width) {
    .w32 => i32,
    .w64 => i64,
};

pub const HUInt = switch (config.highs_int_width) {
    .w32 => u32,
    .w64 => u64,
};

test "public Highs integer types match configured width" {
    switch (config.highs_int_width) {
        .w32 => {
            try std.testing.expect(HInt == i32);
            try std.testing.expect(HUInt == u32);
        },
        .w64 => {
            try std.testing.expect(HInt == i64);
            try std.testing.expect(HUInt == u64);
        },
    }
}
