# !/usr/bin/env ruby
require "test/unit/assertions"
include Test::Unit::Assertions
require "tempfile"

=begin
    Test lib/file_helper.rb:

        - fatal non-existing file

=end

# require .rb files in library
Dir[File.join(__dir__, "..", "lib", "*.rb")].each do |file|
    require file
end

# FileHelper.file_exist_or_die
# - existing file (this file)
assert_nothing_raised {FileHelper.file_exist_or_die(__FILE__)}

# - non-existing file
$stderr.reopen(File.new('/dev/null', 'w')) # suppress SystemExist message
assert_raises SystemExit, FileHelper.file_exist_or_die("foo")
