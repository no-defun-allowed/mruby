MRuby::CrossBuild.new('minimal_plus') do |conf|
  conf.toolchain :gcc
  conf.cc do |cc|
    # Remove some stuff we don't need
    cc.flags = "-Os -fno-exceptions -fno-asynchronous-unwind-tables"
    cc.defines << 'MRB_SMALLER'
  end

  conf.gem 'mrbgems/mruby-bin-mirb/'

  conf.gem 'mrbgems/mruby-bigint/'
  conf.gem 'mrbgems/mruby-compiler/'
  conf.gem 'mrbgems/mruby-io/'
end
