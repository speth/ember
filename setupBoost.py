"""
Compile the necessary Boost libraries for Windows and copy the resulting
files to a location where we can find them.

Run from the main project directory using the version of Python you want
to generate the Boost Python library for.

The script takes as an argument the path to the Boost directory
"""

import sys
import os
import shutil
import platform

from distutils.sysconfig import get_config_var
from os.path import join

def main(boost_dir):
    print 'This is Python %s' % platform.python_version()
    print 'Compiled with %s' % platform.python_compiler()
    print 'Architecture: %s' % platform.architecture()[0]
    buildlibs = ['python', 'filesystem', 'system'] # 'thread', 'date_time'

    pycomp = platform.python_compiler()
    if pycomp.startswith('MSC v.1400'):
        msvc_version = '8.0' # Visual Studio 2005
    elif pycomp.startswith('MSC v.1500'):
        msvc_version  = '9.0' # Visual Studio 2008
    elif pycomp.startswith('MSC v.1600'):
        msvc_version  = '10.0' # Visual Studio 2010
    else:
        print 'Error: Python compiled with unknown version of Visual Studio'
        sys.exit(0)
    msvc_short = msvc_version.replace('.', '')

    assert os.path.exists(boost_dir)
    old_path = os.getcwd()

    python_version = '.'.join(platform.python_version_tuple()[:2])
    python_path = get_config_var('prefix').replace('\\', '/').replace(' ', '\\ ')
    jamConfigText = "using msvc : %s ;\nusing python : %s : %s ;\n" % (msvc_version, python_version, python_path)

    os.chdir(boost_dir)
    jamConfig = open('tmp-config.jam', 'w')
    jamConfig.write(jamConfigText)
    jamConfig.close()

    if not os.path.exists('bjam.exe'):
        os.system('bootstrap.bat')

    boost_version_string = getBoostVersion(boost_dir)
    print 'Boost version: %s' % boost_version_string

    extra_args = ['-j8',
                  '--user-config=%s' % jamConfig.name]
    assert os.path.exists('bjam.exe')

    if '64 bit' in platform.python_compiler():
        extra_args.append('address-model=64')

    extra_args = ' '.join(extra_args)
    libstr = ' '.join('--with-%s' % lib for lib in buildlibs)
    stagedir = join(old_path)

    os.system('bjam.exe %s stage --stagedir=%s runtime-debugging=off variant=release link=shared %s' % (libstr, stagedir, extra_args))
    os.system('bjam.exe %s stage --stagedir=%s runtime-debugging=on variant=debug link=shared %s' % (libstr, stagedir, extra_args))

    os.remove(jamConfig.name)

    os.chdir(old_path)
    for lib in buildlibs:
        for variant,code in (('debug','mt-gd'),
                             ('release','mt')):
            if '64 bit' in platform.python_compiler():
                srcdir = join(boost_dir, 'bin.v2', 'libs', lib, 'build',
                              'msvc-%s' % msvc_version,
                              variant, 'address-model-64', 'threading-multi')
            else:
                srcdir = join(boost_dir, 'bin.v2', 'libs', lib, 'build',
                              'msvc-%s' % msvc_version,
                              variant, 'threading-multi')
            substitutions = (lib, msvc_short, code, boost_version_string)
            libname = 'boost_%s-vc%s-%s-%s.lib' % substitutions
            dllname = 'boost_%s-vc%s-%s-%s.dll' % substitutions
    #        shutil.copy(join(srcdir, libname), join('lib', libname))
    #        print 'successfully copied %s' % libname
    #        shutil.copy(join(srcdir, dllname), join('lib', dllname))
    #        print 'successfully copied %s' % dllname


def getBoostVersion(boost_dir):
    for line in open(join(boost_dir, 'boost', 'version.hpp')):
        if line.startswith('#define BOOST_LIB_VERSION'):
            return line.split('"')[-2]



if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        print r'Usage: setupBoost.py c:\path\to\boost'
