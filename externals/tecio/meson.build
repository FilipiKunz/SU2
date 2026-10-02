check_dir = run_command(python, 
                        script_path / 'check_dir.py',
                        'boost')
if check_dir.returncode() != 0
  message('Extracting boost ...')
  extract_boost = run_command(python, 
                              script_path / 'extract_file.py',
                              'boost.tar.gz', 
                              meson.current_source_dir(), check: true)
else
  message('Boost sources found.')
endif

if mpi
  subdir('teciompisrc')
  su2_deps += teciompi_dep
  if host_machine.system() == 'linux'
    # TecIO's MPI build corrupts large unpartitioned FE zones. Keep its
    # partitioned 3D writer and isolate the serial writer in a shared library.
    subdir('teciosrc')
    tecio_serial = shared_library('su2tecioserial', 'serial_writer.cpp',
      link_with: tecio,
      link_args: ['-Wl,-Bsymbolic', '-Wl,--exclude-libs,ALL'],
      install: true)
    su2_deps += declare_dependency(link_with: tecio_serial)
    su2_cpp_args += '-DHAVE_SERIAL_TECIO'
    su2_install_rpath = '$ORIGIN/..' / get_option('libdir')
  endif
  su2_cpp_args += ['-DHAVE_TECIO', '-DTECIOMPI']
else
  subdir('teciosrc')
  pthread_dep = dependency('threads')
  su2_deps += [tecio_dep, pthread_dep]
  su2_cpp_args += '-DHAVE_TECIO'
endif
