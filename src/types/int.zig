const std = @import("std");
const config = @import("config");

pub const HInt = switch (config.highs_int_width) {
    ._32 => i32,
    ._64 => i64,
};

pub const HUInt = switch (config.highs_int_width) {
    ._32 => u32,
    ._64 => u64,
};

test "public Highs integer types match configured width" {
    switch (config.highs_int_width) {
        ._32 => {
            try std.testing.expect(HInt == i32);
            try std.testing.expect(HUInt == u32);
        },
        ._64 => {
            try std.testing.expect(HInt == i64);
            try std.testing.expect(HUInt == u64);
        },
    }
}
