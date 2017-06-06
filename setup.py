
import sys
from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ["-std=c++11"]

if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS = ["-stdlib=libc++"]

severity = cythonize([
    Extension("severity.weights",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=["severity/weights.pyx",
            "src/weighted_choice.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("severity.simulation",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=["severity/simulation.pyx",
            "src/simulate.cpp",
            "src/weighted_choice.cpp"],
        include_dirs=["src/"],
        language="c++"),
    ])

setup (name="severity",
        description="Package to examine severity of de novo mutations.",
        version="1.0.0",
        author="Jeremy McRae",
        author_email="jeremy.mcrae@sanger.ac.uk",
        license="MIT",
        packages=["severity"],
        install_requires=["denovonear >= 0.4.1",
        ],
        classifiers=[
            "Development Status :: 3 - Alpha",
            "License :: OSI Approved :: MIT License",
        ],
        ext_modules=severity,
        test_suite="tests")
