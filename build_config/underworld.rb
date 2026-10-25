MRuby::CrossBuild.new('underworld') do |conf|
  conf.toolchain :gcc
  DEVKITPRO_PATH = "/opt/devkitpro"
  BIN_PATH = "#{DEVKITPRO_PATH}/devkitPPC/bin"
  conf.gem core: 'hal-wii-io'

  # C compiler
  conf.cc do |cc|
    cc.command = "#{BIN_PATH}/powerpc-eabi-gcc"
    cc.compile_options = %(%{flags} -o "%{outfile}" -c "%{infile}")
    cc.flags = "-Oz -fno-exceptions -fno-asynchronous-unwind-tables"
    cc.defines << 'MRB_SMALLER'
    cc.defines << 'MRB_TARGET_WII'
  end

  # C++ compiler
  conf.cxx do |cxx|
    cxx.command = "#{BIN_PATH}/powerpc-eabi-g++"
    cxx.include_paths = conf.cc.include_paths.dup
    cxx.flags = conf.cc.flags.dup
    cxx.defines = conf.cc.defines.dup
    cxx.compile_options = conf.cc.compile_options.dup
  end

  # Linker
  conf.linker do |linker|
    linker.command = "#{BIN_PATH}/powerpc-eabi-gcc"
  end

  # No executables
  conf.bins = []

  # Do not build executable test
  conf.build_mrbtest_lib_only

  conf.gem 'mrbgems/mruby-bigint/'
  conf.gem 'mrbgems/mruby-compiler/'
  conf.gem 'mrbgems/mruby-enumerator/'
  conf.gem 'mrbgems/mruby-fiber/'
  conf.gem 'mrbgems/mruby-io/'
  conf.gem 'mrbgems/mruby-method/'
  conf.gem 'mrbgems/mruby-pack/'
end
