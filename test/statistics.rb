# !/usr/bin/env ruby
require "test/unit/assertions"
include Test::Unit::Assertions
require "tempfile"

=begin
    Test lib/statistics:

        - mean
        - median
        - fraction
        - percentage

=end

# require .rb files in library (including all subfolders)
Dir[File.join(__dir__, "..", "lib", "**", "*.rb")].each do |file|
    require File.absolute_path(file)
end

# Statistics.mean
# - arithmetic mean of [4, 36, 45, 50, 75] is 42
assert_equal 42, Statistics.mean([4, 36, 45, 50, 75])

# Statistics.median
# - median of [4, 36, 45, 50, 75] is 45
assert_equal 45, Statistics.median([4, 36, 45, 50, 75])
# - median of [4, 36, 39, 45, 50, 75] is 42
assert_equal 42, Statistics.median([4, 36, 39, 45, 50, 75])

# Statistics.fraction
# - float representation of 1/2 is 0.5
assert_equal 0.5, Statistics.fraction(1, 2)

# Statistics.percentage
# - percentage representation of 1/2 is 50
assert_equal 50, Statistics.percentage(1, 2)
