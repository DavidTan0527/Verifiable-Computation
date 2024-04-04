
# Installation Requirements
`Sage >= 10.2` ([Download from Github releases](https://github.com/sagemath/sage/releases))


Ensure that `sage` is added to your `PATH` environment variable by typing `sage -h` in your terminal. You should see something similar to the following prompt:

```
$ sage -h
SageMath version 10.2, Release Date: 2023-12-03

Running Sage:

  file.[sage|py|spyx] -- run given .sage, .py or .spyx file
  -c <cmd>            -- evaluate cmd as sage code
...
```

# Getting Started

## Running the demo code
1. Open a sage session with `sage` at the root directory of the project.
2. Run: `load("load_all.sage"); load("demo.sage");`

## Building your own layer
Refer to `demo.sage` for an example.

# TODO
- [x] Complete utility classes to build the model.
- [x] Complete an all-in-one script.
- [x] Update README with steps.
- [ ] Add documentation for interface methods and usages.
