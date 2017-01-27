# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages


setup(
    name="q2-state-unifrac",
    version="2017.1.0.dev0",
    packages=find_packages(),
    install_requires=['qiime2 == 2017.2.*',
                      'q2-types == 2017.2.*'],
    author="Daniel McDonald",
    author_email="wasade@gmail.com",
    url="https://github.com/wasade/q2-state-unifrac",
    description="An implementation of Strided State UniFrac",
    entry_points={
        "qiime2.plugins":
        ["q2-state-unifrac=q2_state_unifrac.plugin_setup:plugin"]
    },
)
