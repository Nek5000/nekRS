# Contributor's Guidelines to ADIOS 2

This guide will walk you through how to submit changes to ADIOS 2 and interact
with the project as a developer. Information found on ADIOS 2 wiki: https://github.com/ornladios/ADIOS2/wiki under the Contributing to ADIOS section. 

Table of Contents
=================

   * [Contributor's Guide](#contributors-guide)
   * [Table of Contents](#table-of-contents)
      * [Workflow](#workflow)
      * [Setup](#setup)
      * [Making a change and submitting a pull request](#making-a-change-and-submitting-a-pull-request)
         * [Create the topic branch](#create-the-topic-branch)
            * [Do I need to merge master into my branch first?](#do-i-need-to-merge-master-into-my-branch-first)
         * [Submit a pull request](#submit-a-pull-request)
      * [Template implementation separation](#template-implementation-separation)
         * [Example](#example)
            * [Before separation of public and private template implementation](#before-separation-of-public-and-private-template-implementation)
               * [Foo.h](#fooh-containing-all-implementation)
            * [After separation of public and private template implementation](#after-separation-of-public-and-private-template-implementation)
               * [Foo.h](#fooh-containing-only-prototypes-and-explicit-instantiation-declarations)
               * [Foo.inl](#fooinl-containing-template-implementations-that-always-need-to-be-included)
               * [Foo.tcc](#footcc-containing-template-implementations-that-should-be-restricted-to-only-known-types)
               * [Foo.cpp](#foocpp-containing-non-template-implementations-and-explicit-instantiations-definitions-for-known-types)
      * [Code formatting and style](#code-formatting-and-style)

## Workflow
ADIOS uses the GitHub fork-and-branch model. In this, the project "lives" in it's main repository located at https://github.com/ornladios/adios2.git, while each individual developer has their own copy of the repo to work in.  Changes are then submitted to the main repository via pull-requests made with branches from your fork.

## Setup
To setup your local repository for development:

  1. Fork the main repository on GitHub:
     1. Navigate to https://github.com/ornladios/adios2 in your browser.
     1. Click the `[Fork]` button in the upper right-hand side of the page.
  2. Clone the upstream repository to your local machine:
```
$ mkdir adios
$ cd adios
$ git clone https://github.com/ornladios/adios2.git source
Cloning into 'source'...
remote: Counting objects: 4632, done.
remote: Compressing objects: 100% (80/80), done.
remote: Total 4632 (delta 33), reused 0 (delta 0), pack-reused 4549
Receiving objects: 100% (4632/4632), 1.23 MiB | 224.00 KiB/s, done.
Resolving deltas: 100% (2738/2738), done.
Checking connectivity... done.
$
```
  3. Run the `scripts/development/setup.sh` script.  The script will configure an `upstream` remote and link your local master branch to the upstream.
```
$ cd source/
$ ./scripts/developer/setup.sh 
Enter your GitHub username: chuckatkins
Setup SSH push access? [(y)/n] y
Re-configuring local master branch to use upstream
Fetching origin
remote: Counting objects: 6, done.
remote: Compressing objects: 100% (6/6), done.
remote: Total 6 (delta 0), reused 0 (delta 0), pack-reused 0
Unpacking objects: 100% (6/6), done.
From https://github.com/chuckatkins/adios2
Fetching upstream
From https://github.com/ornladios/adios2
 * [new branch]      master     -> upstream/master
 * [new branch]      dashboard  -> upstream/dashboard
 * [new branch]      hooks      -> upstream/hooks
Setting up git aliases...
Setting up git hooks...
$
```

## Making a change and submitting a pull request
At this point you are ready to get to work.  The first thing to do is to create a branch.  ADIOS uses a "branchy" workflow where all changes are committed through self-contained "topic branches".  This helps ensure a clean traceable git history and reduce conflicts.

### Create the topic branch

1. Make sure you are starting from a current master:
```
$ git checkout master
$ git pull
```
2. Create a branch for your change:
```
$ git checkout -b <your-topic-branch-name>
```
3. Make your changes and commits to the branch.
4. Push the branch to your fork:
```
$ git push -u origin HEAD
Counting objects: 189, done.
Delta compression using up to 8 threads.
Compressing objects: 100% (134/134), done.
Writing objects: 100% (189/189), 70.30 KiB | 0 bytes/s, done.
Total 189 (delta 128), reused 102 (delta 44)
remote: Resolving deltas: 100% (128/128), completed with 82 local objects.
To git@github.com:<your-GitHub-username-here>/adios2.git
 * [new branch]      HEAD -> <your-topic-branch-name>
Branch <your-topic-branch-name> set up to track remote branch <your-topic-branch-name> from origin.
$
```

#### Do I need to merge master into my branch first?
Not usually.  The only time to do that is to resolve conflicts.  You're pull request will be automatically rejected if merge-conflicts exist, in which case you can then resolve them by either re-basing your branch the current master (preferable):
```
$ git fetch --all -p
...
$ git rebase upstream/master
$ git push -f
```
or if necessary or re-basing is not a viable option then you can always fall back to merging in master but it should be generally discouraged as it tends to make the git history difficult to follow:
```
$ git fetch --all -p
...
$ git merge upstream/master
$ git push -f
```

### Submit a pull request
1. Log in to your GitHub fork.
2. You should see a message at the top that informs you of your recently pushed branch, something like: `<your-topic-branch-name> (2 minutes ago)`.  On the right side, select the `[Compare & pull request]` button.
3. Fill in the appropriate information for the name of the branch and a brief summary of the changes it contains.
   * The default configuration will be for the topic branch to be merged into the upstream's master branch.  You can change this if you want to submit to a different branch.
4. Click `[Create pull request]`.

You have now created a pull request (PR) that is pending several status checks before it can be merged.  Currently, the only check being performed is for source code formatting and style.  In the future, however, the will be a more in depth continuous integration system tied to the pull requests that tests for build and test failures every time a PR is submitted or updated.  Once the status checks pass, the PR will be eligible for merging by one of the project maintainers.

## Template implementation separation
The ADIOS C++ classes try to explicitly separate class declarations from their implementation.  Typically this is done by having a separate .h and .cpp file,  however it get's more complicated when templates are involved.  To maintain the distinct separation between definition and implementation, we use explicit instantiation with 4 different source file types:
* ClassName.h
  * The main header file containing *only* the class and member declarations with no implementation.  This also contains the declarations for explicitly instantiated members.
* ClassName.inl
  * A file containing inline function implementations that need to be made public.  This is to be included at the bottom of ClassName.h and should *only* contain implementations that need to be made public.
* ClassName.tcc
  * A file containing most of the template implementations that can be hidden through explicit instantiation.
* ClassName.cpp
  * A file containing the non-template implementations and the explicit instation of any template members.

### Example
Here is an example of a simple class `Foo` with template member functions `Bar1` and `Bar2`

#### Before separation of public and private template implementation
##### Foo.h containing all implementation
```cpp
#ifndef FOO_H_
#define FOO_H_

namespace adios
{

class Foo
{
public:
    Foo()
    : m_Bar1Calls(0), m_Bar2Calls(0), m_Bar3Calls(0);
    {
    }

    virtual ~Foo() = default;

    template<typename T>
    void Bar1()
    {
        Bar1Helper<T>();
    }

    template<typename T>
    void Bar2()
    {
        Bar2Helper<T>();
    }

    void Bar3()
    {
        Bar3Helper();
    }

private:
    template<typename T>
    void Bar1Helper()
    {
        ++m_Bar1Calls;
    }

    template<typename T>
    void Bar2Helper()
    {
        ++m_Bar2Calls;
    }

    void Bar3Helper()
    {
        ++m_Bar3Calls;
    }

    size_t m_Bar1Calls;
    size_t m_Bar2Calls;
    size_t m_Bar3Calls;
};

} // end namespace adios
#endif // FOO_H_
```

#### After separation of public and private template implementation

In this example, we want to hide the template implementation from the header.  We will implement this such that `Bar1` is only callable from the core numeric types, i.e. ints, floats, and complex, while `Bar2` is callable from all types.  This will necessitate that `Bar1` and it's helper function is implemented in a .tcc file with explicit instantiation for the allowed types while `Bar2` and it's helper function will need to be inlined in the .inl file to be accessible for all types.  We will also use a helper macro ADIOS provides to iterate over the core numeric types for the explicit instantiation of `Bar1`.

##### Foo.h containing only prototypes and explicit instantiation declarations
```cpp
#ifndef FOO_H_
#define FOO_H_

#include "ADIOSMacros.h"

namespace adios
{
class Foo
{
public:
    Foo();
    virtual ~Foo() = default;

    template<typename T>
    void Bar1();

    template<typename T>
    void Bar2();

    void Bar3();
private:
    template<typename T>
    void Bar1Helper();

    template<typename T>
    void Bar2Helper();

    void Bar3Helper;

    size_t m_Bar1Calls;
    size_t m_Bar2Calls;
    size_t m_Bar3Calls;
};

// Create declarations for explicit instantiations
#define declare_explicit_instantiation(T)       \
    extern template void Foo::Bar1<T>();

ADIOS_FOREACH_STDTYPE_1ARG(declare_explicit_instantiation)
#undef(declare_explicit_instantiation)
} // end namespace adios

#include "Foo.inl"
#endif // FOO_H_
```
Note here that Bar1Helper does not need an explicit instantiation because it's not a visible funtion in the callable interface.  It's implementaion will be available to Bar1 inside the tcc file where it's called from.

##### Foo.inl containing template implementations that always need to be included
```cpp
#ifndef FOO_INL_
#define FOO_INL_
#ifndef FOO_H_
#error "Inline file should only be included from it's header, never on it's own"
#endif

// No need to include Foo.h since it's where this is include from

namespace adios
{

template<typename T>
void Foo::Bar2()
{
    Bar2Helper<T>();
}

template<typename T>
void Foo::Bar2Helper()
{
    ++m_Bar2Calls;
}

} // end namespace adios

#endif // FOO_INL_
```

##### Foo.tcc containing template implementations that should be restricted to only known types
```cpp
#ifndef FOO_TCC_
#define FOO_TCC_

#include "Foo.h"
namespace adios
{

template<typename T>
void Foo::Bar1()
{
    Bar1Helper<T>();
}

template<typename T>
void Foo::Bar1Helper()
{
    ++m_Bar1Calls;
}

} // end namespace adios

#endif // FOO_TCC_
```

##### Foo.cpp containing non-template implementations and explicit instantiations definitions for known types.
```cpp
#include "Foo.h"
#include "Foo.tcc"

namespace adios
{

Foo::Foo()
: m_Bar1Calls(0), m_Bar2Calls(0), m_Bar3Calls(0)
{
}

void Foo::Bar3()
{
    Bar3Helper();
}

void Foo::Bar3Helper()
{
    ++m_Bar3Calls;
}

// Create explicit instantiations of existing definitions
#define define_explicit_instantiation(T)  \
    template void Foo::Bar1<T>();

ADIOS_FOREACH_STDTYPE_1ARG(define_explicit_instantiation)
#undef(define_explicit_instantiation)

} // end namespace adios
```

## Code formatting and style using clang-format 7
ADIOS uses the [clang-format version 7](https://releases.llvm.org/7.0.0/tools/clang/docs/ClangFormat.html) tool to automatically enforce source code style and formatting rules.  There are various ways to integrate the clang-format tool into your IDE / Code Editor depending on if you use Emacs, Vim, Eclipse, KDevelop, Microsoft Visual Studio, etc. that are a bit outside the scope of this document. A quick google search for "integrate <insert-editor-here> clang-format" should point you in the right direction.  However, you can always reformat the code manually by running:

```
clang-format -i SourceFile.cpp SourceFile.h
```

That will apply the formatting rules used by the ADIOS project.

While some of the formatting rules are fairly detailed, the main points are:

1. Lines no longer than 80 characters.
1. Always use braces { and }, even for 1 line if blocks.
1. Use 4 spaces for indentation.

There are more formatting rules but these three should at least get you close and prevent any drastic re-writes from the re-formatting tools.  More details can be found by looking at the `.clang-format` config file in the root of the repository and by looking at the [clang-format documentation](https://releases.llvm.org/7.0.0/tools/clang/docs/ClangFormat.html).
