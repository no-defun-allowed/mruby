MRuby::Gem::Specification.new('hal-wii-io') do |spec|
  spec.license = 'MIT'
  spec.author  = 'Hayley Patton'
  spec.summary = 'POSIX HAL for libogc'

  # HAL gem depends on feature gem - brings in mruby-io automatically
  spec.add_dependency 'mruby-io', core: 'mruby-io'
end
