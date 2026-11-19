#!/usr/bin/env ruby

# Usage: ruby install.rb [clean] [name=value ...]
#
# Builds the spherical-coordinate 3D executables x_eikonal3, x_fatt3 and
# x_tloc3 in the foreground.
#
# clean removes the old objects and executables before compiling. Each
# name=value goes to make and overrides the variable of the same name in the
# Makefiles, e.g.,
#
#   ruby install.rb clean flitdir=$HOME/flit

require 'shellwords'

clean = ARGV.include?("clean")
options = ARGV.grep(/\A\w+=/)
unknown = ARGV - ["clean"] - options
if not unknown.empty?
        abort("Unknown argument: #{unknown.join(' ')}\nUsage: ruby install.rb [clean] [name=value ...]")
end
opts = Shellwords.join(options)

# The spherical branch has 3D executables only
if clean
        system(" cd 3d; make clean #{opts}; make #{opts}")
else
        system(" cd 3d; make #{opts}")
end
