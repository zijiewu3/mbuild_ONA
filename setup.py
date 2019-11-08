from setuptools import setup

setup(
      name="mbuild_ONA",
      install_requires="mbuild",
      entry_points={
                    "mbuild.plugins":[ "ONA_box = mbuild_ONA.mbuild_ONA:ONA_box"
                        ]
                    },
                    py_modules=["mbuild_ONA"],
                        )
