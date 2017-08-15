# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

import subprocess
import os


SUCPP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'sucpp/')


def compile_ssu():
    """Clean and compile the SSU binary"""
    # clean the target
    subprocess.call(['make', 'clean'], cwd=SUCPP)

    cmd = ['make', 'test']
    ret = subprocess.call(cmd, cwd=SUCPP)
    if ret != 0:
        raise Exception('Error compiling ssu!')

    cmd = ['make', 'main']
    ret = subprocess.call(cmd, cwd=SUCPP)
    if ret != 0:
        raise Exception('Error compiling ssu!')


class PostBuildCommand(install):
    """Post-installation for development mode."""
    def run(self):
        install.run(self)
        self.execute(compile_ssu, [], 'Compiling SSU')
        self.copy_file(os.path.join(SUCPP, 'ssu'),
                       os.path.join(self.install_libbase, 'q2_state_unifrac/'))


class PostDevelopCommand(develop):
    """Post-installation for development mode (i.e. `pip install -e ...`)."""
    def run(self):
        develop.run(self)
        self.execute(compile_ssu, [], 'Compiling SSU')
        self.copy_file(os.path.join(SUCPP, 'ssu'),
                       os.path.join(self.install_libbase, 'q2_state_unifrac/'))


setup(
    name="q2-state-unifrac",
    version="2017.2.0",
    packages=find_packages(),
    author="Daniel McDonald",
    license='BSD-3-Clause',
    author_email="wasade@gmail.com",
    url="https://github.com/wasade/q2-state-unifrac",
    description="An implementation of Strided State UniFrac",
    entry_points={
        "qiime2.plugins":
        ["q2-state-unifrac=q2_state_unifrac.plugin_setup:plugin"]
    },
    cmdclass={'install': PostBuildCommand, 'develop': PostDevelopCommand}
)
