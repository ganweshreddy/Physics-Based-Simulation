# Style and Documentation Guidelines

Our coding style should be as clean and easy to read as possible.
To assist with style guide compliance a `.clang-format` file is included in this repository.

## General

### C++ Version

Currently, code should target C++17. The targeted version will be changed over time.

### Files

Every C/C++ implementation file should have the extension `*.cpp`. Besides some common exceptions they should always have an associated header file with extension `*.hpp`.
Similar rules apply to CUDA files with extensions `*.cu` and `*.cuh` respectively.

Generally, every class should have its own file. The filename should always reflect that. Filenames have to be snake case. Header files only containing forward declarations should add the suffix `'_fwd'`.

The following example illustrates these conventions:

`base_array.hpp`
```C++
namespace vislab
{
    class BaseArray
    {
    public:
        BaseArray();
    };
}
```

`base_array.cpp`
```C++
using namespace vislab;

BaseArray::BaseArray()
{
}
```

`base_array_fwd.hpp`
```C++
namespace vislab
{
    class BaseArray;
}
```

### Headers

All header files should be guarded to prevent multiple inclusion. The `#pragma once` directive must always be used to guard the header file. The classic guard defines can be provided *in addition* to that but are not required.

Include only what you need, while refusing to rely on transitive inclusions.

### Namespaces

Place code in namespaces. The common namespace in this project is the `vislab` namespace. Avoid *using-directives* in header files.

### Classes

Define all non-trivial constructors and operators when they are useful and relevant.
Data members should generally kept private. Access can be granted by public get/set member functions.

Similar declarations should be grouped together.
Usually, any class definition starts with a `public:` section followed by the `protected:` and then the `private:` sections. Empty sections should be omitted.

Short member functions may be defined inline. Larger definitions should be placed outside the class declaration.

### Naming

Names should *always* be reasonably short and descriptive!

## Syntax

Generally a camel case convention is used.

### Classes

Classes, structs and other types are upper camel case.

```C++
class VectorField;
using Matrix3d = MatrixBase<double, 3, 3>;
```

Additional prefixes should be introduced for enumerations (``E``) and interfaces (``I``).

```C++
class IField; // Interface for fields.
enum class EFieldType; // Enumeration specifying field types.
```

Private member variables should be prefixed with ``m``:

```C++
int mIndex;
double mAngularVelocity;
```

Public member variables should generally follow the style of functions.

### Functions

Functions as well as public members use the lower camel case or dromedary case convention.

```C++
double getAngularVelocity();
void setIndex(int idx);
```

### Variables

Local Variables also use lower came case.

```C++
double volumetricDensity;
std::string namedAttribute;
```

### Additional

Everything else is done automatically by using clang-format and the configuration file included in the repository.

This includes among some other things:
* Space indentation
* Breaks before braces
* Breaks in inheritance lists and constructor initializers
* Spaces in braced lists
* Alignment of consecutive assignments

Formatting has to be applied at least before pushing to the repository to guarantee some consistency.

To enable clang-format in visual studio see the reference page [here](https://learn.microsoft.com/en-us/visualstudio/ide/reference/options-text-editor-c-cpp-formatting?view=vs-2022).
In visual studio code the setup is easy as well. By default the [C/C++ Extension by Microsoft](https://code.visualstudio.com/docs/languages/cpp)  is able to [run](https://code.visualstudio.com/docs/cpp/cpp-ide) clang-format. Other extensions can be used as well. Other IDE's or editors are usually able to integrate clang-format as well; nevertheless command line is still an option if there is a strong aversion to anything with a GUI.

## Documentation

Documentation helps others and yourself better understand your code, and is therefore essential.
Doxygen is used to generate readable API documentation automatically from source code comments. Writing those comments wherever possible is important and will be reviewed.

Doxygen comments should be written in the Default Doxygen JavaDoc-style (C-style) since this style provides the best general support and can be seen as the de-facto standard.

Function documentation is probably the most important part here. A brief description with `@brief` is mandatory. For longer descriptions or important details `@details` is useful as well. Template parameters should be described by `@tparam`. Function parameters must be documented with `@param` - input and output parameters should be documented as well. Documentation of the return value with `@return` is mandatory. The following example showcases the baseline.

```C++
/**
 * @brief The general description of the function.
 * @details Important details to mention.
 * @tparam T A template parameter with symbol T.
 * @param[in] value An input parameter with the name 'value'.
 * @param[out] executed An output parameter with the name 'executed'.
 * @param[in,out] bitMask An input and output parameter with the name 'bitMask'.
 * @return The return parameter.
 */
template <typename T>
int anyFunction(double value, bool executed, int bitMask);

/** @brief Member variable description.*/
double mVariable;
```


Inline documentation can be used as well but should be avoided for better readability.

```C++
/**
 * @brief General description of the enumeration.
 */
enum class EType
{
	UNSIGNED, /**< @brief Inline description. */
	RGB,      /**< @brief Inline description. */
	TEAPOT    /**< @brief Inline description. */
};
```

Useful additions to emphasize certain aspects of the documentation:

```C++
/**
 * Italic formatting with @a to reference other instances
 * Bold formatting with @b to emphasize words.
 * Monospace (typewriter) formatting with @c for file names/code symbols.
 */
```

### Documentation Guide

The following documentation is the required baseline to pass the code reviews.

### Files

```C++
/**
 * @file      [Mandatory]
 * @author    [Mandatory]
 * @version   [Mandatory]
 * @date      [Mandatory]
 * @copyright [Mandatory]
 */
```

### Functions

```C++
/**
 * @brief [Mandatory]
 * @details [Optional]
 * @tparam T [Mandatory if existent]
 * @param[in,out] value [Mandatory if existent]
 * @return  [Mandatory if existent]
 */
```

### Variables

```C++
/**
 * @brief [Mandatory]
 * @details [Optional]
 */
```

### Classes

```C++
/**
 * @brief [Mandatory]
 * @details [Optional]
 */
```

### Structs

```C++
/**
 * @brief [Mandatory]
 * @details [Optional]
 */
```

### Enumerations

```C++
/**
 * @brief [Mandatory]
 * @details [Optional]
 */
```

### Unions

```C++
/**
 * @brief [Mandatory]
 * @details [Optional]
 */
```

### Namespaces

```C++
/**
 * @brief [Mandatory]
 * @details [Optional]
 */
```

### Defines / Typedef / Using

```C++
/**
 * @brief [Mandatory]
 * @details [Optional]
 */
```

Usually all declarations in header files (``*.hpp``) should be documented - especially API. In implementation files (``*.cpp``) the file documentation is enough most of the time.
Doxygen by default warns about undocumented parts. This usually is a good guideline to decide, what should be documented and what can be left out. In the rare case this warning does not fit, a part can be left out by:

```C++
//** @cond NO_COND */

int undocumented;

//** @endcond */
```

Keep in mind: This conditionally only disables the parsing of that part of the file but can be enabled nonetheless when checking the documentation in reviews.

### Auto generate some comments

Visual Studio has [Doxygen Comment Support](https://devblogs.microsoft.com/cppblog/doxygen-and-xml-doc-comment-support/) since Visual Studio 2019 16.6.2. The comment style can be setup there.

In visual studio code the [C/C++ Extension by Microsoft](https://code.visualstudio.com/docs/languages/cpp) does already include basic comment generation. If we need to customize more then we can to match the Visual Studio generated comments, we can switch to the [Doxygen Documentation Generator](https://marketplace.visualstudio.com/items?itemName=cschlosser.doxdocgen) extension. The latter is a lot more customizable and allows for more settings regarding i.e. file documentation.
Other IDEs and editors usually have options as well, manually writing Doxygen comments is of course possible as well.

## Conclusion

Readable code, understandable documentation, and consistency are extremely important, especially for larger projects with long lifespans. Following the guidelines helps immensely. If you get used to adapting your own standard to this end, it increases the quality of the entire software.