flag_list = Split("""
                  -march=x86-64
                  -pedantic
                  -std=c99 
                  -g
                  """)
env = Environment(SHLIBSUFFIX = '.so', CCFLAGS = flag_list, CPPPATH = '../include')
Export('env')

nsatmo = SConscript('src/SConscript', variant_dir='./build', duplicate=0)
Default(nsatmo)


env.Install('./lib', nsatmo)
env.Alias('install', ['./lib'])
