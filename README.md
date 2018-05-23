<p align="center">
  <h1 align="center">
    pd.dummies
  </h1>
  <p align="center">
    A set of dummy pure-data objects
  </p>
  <p align="center">
    <a href="https://travis-ci.org/pierreguillot/pd.dummies"><img src="https://img.shields.io/travis/pierreguillot/pd.dummies.svg?label=travis" alt="Travis CI"></a>
    <a href="https://ci.appveyor.com/project/pierreguillot/pd.dummies/history"><img src="https://img.shields.io/appveyor/ci/pierreguillot/pd-dummies.svg?label=appveyor" alt="Appveyor CI"></a>
  </p>
</p>

### Presentation

The dummies objects are a set of examples used for teaching purposes or to help other developers. You should not rely on it for your creations or any project that aims to be sustainable. By the way, most of the objects should be pretty useless.

### Dependencies

- [pure-data](https://github.com/pure-data/pure-data.git)
- [pd.build](https://github.com/pierreguillot/pd.build.git)

## Compilation

```
git submodule update --init --recursive
mkdir build && cd build
cmake ..
cmake --build . --config Release
```
