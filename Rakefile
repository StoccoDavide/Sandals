# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Sandals project is distributed under the GNU GPLv3.                     #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# e-mail: davide.stocco@unitn.it           e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

# Configuration of the build
BUILD_DEBUG      = false
BUILD_TESTS      = true
BUILD_EXAMPLES   = false
BUILD_BENCHMARKS = false

case RUBY_PLATFORM
when /mingw|mswin/
  PARALLEL = ''
  QUIET    = ''
else
  require 'etc'
  cmakeversion = %x( cmake --version ).scan(/\d+\.\d+\.\d+/).last
  mm = cmakeversion.split('.');
  if mm[0].to_i > 3 || (mm[0].to_i == 3 && mm[1].to_i >= 12) then
    PARALLEL = "--parallel #{Etc.nprocessors} "
    QUIET    = '-- --quiet '
  else
    PARALLEL = ''
    QUIET    = ''
  end
end

if (/cygwin|mswin|mingw|bccwin|wince|emx/ =~ RUBY_PLATFORM) != nil then
  # Linux
  task :default => [:build_linux]
elsif (/darwin/ =~ RUBY_PLATFORM) != nil then
  # OsX
  task :default => [:build_osx]
else
  # Windows
  task :default => [:build_windows]
end

file_base = File.expand_path(File.dirname(__FILE__)).to_s

cmd_cmake_build = "-G Ninja"
if BUILD_TESTS then
  cmd_cmake_build += ' -DBUILD_TESTS:VAR=true '
else
  cmd_cmake_build += ' -DBUILD_TESTS:VAR=false '
end
if BUILD_EXAMPLES then
  cmd_cmake_build += ' -DBUILD_EXAMPLES:VAR=true '
else
  cmd_cmake_build += ' -DBUILD_EXAMPLES:VAR=false '
end
if BUILD_BENCHMARKS then
  cmd_cmake_build += ' -DBUILD_BENCHMARKS:VAR=true '
else
  cmd_cmake_build += ' -DBUILD_BENCHMARKS:VAR=false '
end
if BUILD_DEBUG then
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=STATUS '
else
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=STATUS '
end

task :default => [:build]

TESTS = []

desc "run tests"
task :run do
  puts "run test".yellow
  Dir.glob('bin/*') do |cmd|
    next if cmd =~ /.manifest$|.dSYM$/
    puts "execute: #{cmd}".yellow
    sh cmd
  end
end

task :build_gen do

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = "cmake " + cmd_cmake_build

  puts "run CMake for Sandals".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMake for Sandals".yellow
  if BUILD_DEBUG then
    sh 'cmake --build . --config Debug ' + PARALLEL + QUIET
  else
    sh 'cmake --build . --config Release ' + PARALLEL + QUIET
  end
  FileUtils.cd '..'
end

desc 'compile for OsX'
task :build_osx => :build_gen do |t, args|
end

desc 'compile for Linux'
task :build_linux => :build_gen do |t, args|
end

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  args.with_defaults( :year => "2017", :bits => "x64" )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = win_vs(args.bits,args.year) + cmd_cmake_build

  puts "run CMake for Sandals".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMake for Sandals".yellow
  if BUILD_DEBUG then
    sh 'cmake --build . --config Debug ' + PARALLEL + QUIET
  else
    sh 'cmake  --build . --config Release  ' + PARALLEL + QUIET
  end

  FileUtils.cd '..'
end

desc "build for OsX/Linux/Windows"
task :build do
  if (/cygwin|mswin|mingw|bccwin|wince|emx/ =~ RUBY_PLATFORM) != nil then
    # Linux
    Rake::Task["build_linux"].invoke
  elsif (/darwin/ =~ RUBY_PLATFORM) != nil then
    # OsX
    Rake::Task["build_osx"].invoke
  else
    # Windows
    Rake::Task["build_windows"].invoke
  end
end

task :clean_gen do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
end

desc "clean for OsX"
task :clean_osx => :clean_gen do
end

desc "clean for Linux"
task :clean_linux => :clean_gen do
end

desc "clean for Windows"
task :clean_win do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  FileUtils.rm_rf 'vs_*'
end

desc "clean for OsX/Linux/Windows"
task :clean do
  if (/cygwin|mswin|mingw|bccwin|wince|emx/ =~ RUBY_PLATFORM) != nil then
    # Linux
    Rake::Task["clean_linux"].invoke
  elsif (/darwin/ =~ RUBY_PLATFORM) != nil then
    # OsX
    Rake::Task["clean_linux"].invoke
  else
    # Windows
    Rake::Task["clean_linux"].invoke
  end
end
