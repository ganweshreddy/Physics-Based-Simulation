#include "vislab/core/numeric_parameter.hpp"

#include "gtest/gtest.h"

namespace vislab
{
    TEST(core, parameter)
    {
        // create a parameter
        auto floatParameter = std::make_unique<FloatParameter>();

        // define an event listener listening to parameter changes and record if a change was reported
        float listened = 0;
        floatParameter->onChange += [=, &listened](Parameter<float>* sender, const float* value)
        {
            listened = sender->getValue();
        };

        // set a testing value
        const float testValue = 1.2f;
        floatParameter->setValue(testValue);

        // see if the event fired successfully and whether the getter gets the correct result
        EXPECT_EQ(testValue, listened);
        EXPECT_EQ(testValue, floatParameter->getValue());

        // unhook the listener
        floatParameter->onChange.clearAll();

        // set another test value
        const float testValue2 = 3.2f;
        floatParameter->setValue(testValue2);

        // see if no event was fired and whether the getter gets the correct result
        EXPECT_NE(testValue2, listened);
        EXPECT_EQ(testValue2, floatParameter->getValue());
    }
}
