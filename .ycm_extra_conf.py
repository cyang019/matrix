
def Settings( **kwargs ):
  return {
    'flags': ['-x', 'c++', '-Wall', '-Wextra', '-Werror',
              '-isystem',
              '-std=++17',
              '-Wc++17-extensions',
              '-I./matrix_impl/include',
              '-I./benchmark_matrix/include',
              '-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers',
              '-I/usr/include',
              '-I./build',
              '-I../build',
              '-I/build/googletest/googletest-src/googletest/include',
              '-I../build/googletest/googletest-src/googletest/include'],
  }
