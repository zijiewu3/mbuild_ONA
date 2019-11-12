from setuptools import setup

setup(
    name="mbuild_ona",
    install_requires="mbuild",
    entry_points={"mbuild.plugins": ["DNA = mbuild_ona.mbuild_ona:DNA"]},
    py_modules=["mbuild_ona"],
)
