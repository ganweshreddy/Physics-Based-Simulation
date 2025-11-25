#include "vislab/core/array.hpp"

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    /**
     * @brief Test function for arrays with variable number of components.
     * @tparam TArray Array type to run the tests for.
     */
    template <typename TArray>
    void test_ArrayX()
    {
        auto myArray = std::make_unique<TArray>();
        myArray->setNumberOfComponents(3);
        Eigen::Index numComponents = myArray->getNumberOfComponents();
        EXPECT_EQ(numComponents, 3);

        myArray->setSize(3);
        EXPECT_EQ(myArray->getSize(), 3);

        myArray->setName("test string");
        EXPECT_EQ(myArray->getName(), "test string");

        // regular setter and getter functions
        myArray->setValue(0, Eigen::Vector<typename TArray::Scalar, 3>(1.1f, 2.1f, 3.1f));
        myArray->setValue(1, Eigen::Vector<typename TArray::Scalar, 3>(1.2f, 2.2f, 3.2f));
        myArray->setValue(2, Eigen::Vector<typename TArray::Scalar, 3>(1.3f, 2.3f, 3.3f));
        Eigen::Vector<typename TArray::Scalar, 3> val0 = myArray->getValue(0);
        Eigen::Vector<typename TArray::Scalar, 3> val1 = myArray->getValue(1);
        Eigen::Vector<typename TArray::Scalar, 3> val2 = myArray->getValue(2);
        EXPECT_EQ(val0.x(), (typename TArray::Scalar)1.1f);
        EXPECT_EQ(val0.y(), (typename TArray::Scalar)2.1f);
        EXPECT_EQ(val0.z(), (typename TArray::Scalar)3.1f);
        EXPECT_EQ(val1.x(), (typename TArray::Scalar)1.2f);
        EXPECT_EQ(val1.y(), (typename TArray::Scalar)2.2f);
        EXPECT_EQ(val1.z(), (typename TArray::Scalar)3.2f);
        EXPECT_EQ(val2.x(), (typename TArray::Scalar)1.3f);
        EXPECT_EQ(val2.y(), (typename TArray::Scalar)2.3f);
        EXPECT_EQ(val2.z(), (typename TArray::Scalar)3.3f);

        // double setter and getter functions
        myArray->setValueDouble(0, as_double_t<Eigen::Vector<typename TArray::Scalar, 3>>(11.1, 12.1, 13.1));
        myArray->setValueDouble(1, as_double_t<Eigen::Vector<typename TArray::Scalar, 3>>(11.2, 12.2, 13.2));
        myArray->setValueDouble(2, as_double_t<Eigen::Vector<typename TArray::Scalar, 3>>(11.3, 12.3, 13.3));
        as_double_t<Eigen::Vector<typename TArray::Scalar, 3>> val0d = myArray->getValueDouble(0);
        as_double_t<Eigen::Vector<typename TArray::Scalar, 3>> val1d = myArray->getValueDouble(1);
        as_double_t<Eigen::Vector<typename TArray::Scalar, 3>> val2d = myArray->getValueDouble(2);
        EXPECT_EQ(val0d.x(), (typename TArray::Scalar)11.1);
        EXPECT_EQ(val0d.y(), (typename TArray::Scalar)12.1);
        EXPECT_EQ(val0d.z(), (typename TArray::Scalar)13.1);
        EXPECT_EQ(val1d.x(), (typename TArray::Scalar)11.2);
        EXPECT_EQ(val1d.y(), (typename TArray::Scalar)12.2);
        EXPECT_EQ(val1d.z(), (typename TArray::Scalar)13.2);
        EXPECT_EQ(val2d.x(), (typename TArray::Scalar)11.3);
        EXPECT_EQ(val2d.y(), (typename TArray::Scalar)12.3);
        EXPECT_EQ(val2d.z(), (typename TArray::Scalar)13.3);

        // int64_t setter and getter functions
        myArray->setValueInt64(0, as_int64_t<Eigen::Vector<typename TArray::Scalar, 3>>(2, 3, 4));
        myArray->setValueInt64(1, as_int64_t<Eigen::Vector<typename TArray::Scalar, 3>>(3, 4, 5));
        myArray->setValueInt64(2, as_int64_t<Eigen::Vector<typename TArray::Scalar, 3>>(4, 5, 6));
        as_int64_t<Eigen::Vector<typename TArray::Scalar, 3>> val0i = myArray->getValueInt64(0);
        as_int64_t<Eigen::Vector<typename TArray::Scalar, 3>> val1i = myArray->getValueInt64(1);
        as_int64_t<Eigen::Vector<typename TArray::Scalar, 3>> val2i = myArray->getValueInt64(2);
        EXPECT_EQ(val0i.x(), (typename TArray::Scalar)2);
        EXPECT_EQ(val0i.y(), (typename TArray::Scalar)3);
        EXPECT_EQ(val0i.z(), (typename TArray::Scalar)4);
        EXPECT_EQ(val1i.x(), (typename TArray::Scalar)3);
        EXPECT_EQ(val1i.y(), (typename TArray::Scalar)4);
        EXPECT_EQ(val1i.z(), (typename TArray::Scalar)5);
        EXPECT_EQ(val2i.x(), (typename TArray::Scalar)4);
        EXPECT_EQ(val2i.y(), (typename TArray::Scalar)5);
        EXPECT_EQ(val2i.z(), (typename TArray::Scalar)6);

        // resize
        myArray->setSize(4);
        EXPECT_EQ(myArray->getSize(), 4);
    }

    /**
     * @brief Test function for arrays with one component.
     * @tparam TArray Array type to run the tests for.
     */
    template <typename TArray>
    void test_Array1()
    {
        auto myArray               = std::make_unique<TArray>();
        Eigen::Index numComponents = myArray->getNumberOfComponents();
        EXPECT_EQ(numComponents, 1);

        myArray->setSize(3);
        EXPECT_EQ(myArray->getSize(), 3);

        myArray->setName("test string");
        EXPECT_EQ(myArray->getName(), "test string");

        // regular setter and getter functions
        myArray->setValue(0, typename TArray::Element(1.1f));
        myArray->setValue(1, typename TArray::Element(1.2f));
        myArray->setValue(2, typename TArray::Element(1.3f));
        typename TArray::Element val0 = myArray->getValue(0);
        typename TArray::Element val1 = myArray->getValue(1);
        typename TArray::Element val2 = myArray->getValue(2);
        EXPECT_EQ(val0.x(), (typename TArray::Scalar)1.1f);
        EXPECT_EQ(val1.x(), (typename TArray::Scalar)1.2f);
        EXPECT_EQ(val2.x(), (typename TArray::Scalar)1.3f);

        // special setter function for scalars
        myArray->setValue(0, 11.1f);
        myArray->setValue(1, 11.2f);
        myArray->setValue(2, 11.3f);
        typename TArray::Element vals0 = myArray->getValue(0);
        typename TArray::Element vals1 = myArray->getValue(1);
        typename TArray::Element vals2 = myArray->getValue(2);
        EXPECT_EQ(vals0.x(), (typename TArray::Scalar)11.1f);
        EXPECT_EQ(vals1.x(), (typename TArray::Scalar)11.2f);
        EXPECT_EQ(vals2.x(), (typename TArray::Scalar)11.3f);

        // double setter and getter functions
        myArray->setValueDouble(0, 21.1);
        myArray->setValueDouble(1, 21.2);
        myArray->setValueDouble(2, 21.3);
        myArray->setValueDouble(0, as_double_t<typename TArray::Element>(21.1));
        myArray->setValueDouble(1, as_double_t<typename TArray::Element>(21.2));
        myArray->setValueDouble(2, as_double_t<typename TArray::Element>(21.3));
        as_double_t<typename TArray::Element> val0d = myArray->getValueDouble(0);
        as_double_t<typename TArray::Element> val1d = myArray->getValueDouble(1);
        as_double_t<typename TArray::Element> val2d = myArray->getValueDouble(2);
        EXPECT_EQ(val0d.x(), (typename TArray::Scalar)21.1);
        EXPECT_EQ(val1d.x(), (typename TArray::Scalar)21.2);
        EXPECT_EQ(val2d.x(), (typename TArray::Scalar)21.3);

        // int64_t setter and getter functions
        myArray->setValueInt64(0, 2);
        myArray->setValueInt64(1, 3);
        myArray->setValueInt64(2, 4);
        myArray->setValueInt64(0, as_int64_t<typename TArray::Element>(2));
        myArray->setValueInt64(1, as_int64_t<typename TArray::Element>(3));
        myArray->setValueInt64(2, as_int64_t<typename TArray::Element>(4));
        as_int64_t<typename TArray::Element> val0i = myArray->getValueInt64(0);
        as_int64_t<typename TArray::Element> val1i = myArray->getValueInt64(1);
        as_int64_t<typename TArray::Element> val2i = myArray->getValueInt64(2);
        EXPECT_EQ(val0i.x(), (typename TArray::Scalar)2);
        EXPECT_EQ(val1i.x(), (typename TArray::Scalar)3);
        EXPECT_EQ(val2i.x(), (typename TArray::Scalar)4);

        // resize
        myArray->setSize(4);
        EXPECT_EQ(myArray->getSize(), 4);
    }

    /**
     * @brief Test function for arrays with two components.
     * @tparam TArray Array type to run the tests for.
     */
    template <typename TArray>
    void test_Array2()
    {
        auto myArray               = std::make_unique<TArray>();
        Eigen::Index numComponents = myArray->getNumberOfComponents();
        EXPECT_EQ(numComponents, 2);

        myArray->setSize(3);
        EXPECT_EQ(myArray->getSize(), 3);

        myArray->setName("test string");
        EXPECT_EQ(myArray->getName(), "test string");

        // regular setter and getter functions
        myArray->setValue(0, typename TArray::Element(1.1f, 2.1f));
        myArray->setValue(1, typename TArray::Element(1.2f, 2.2f));
        myArray->setValue(2, typename TArray::Element(1.3f, 2.3f));
        typename TArray::Element val0 = myArray->getValue(0);
        typename TArray::Element val1 = myArray->getValue(1);
        typename TArray::Element val2 = myArray->getValue(2);
        EXPECT_EQ(val0.x(), (typename TArray::Scalar)1.1f);
        EXPECT_EQ(val0.y(), (typename TArray::Scalar)2.1f);
        EXPECT_EQ(val1.x(), (typename TArray::Scalar)1.2f);
        EXPECT_EQ(val1.y(), (typename TArray::Scalar)2.2f);
        EXPECT_EQ(val2.x(), (typename TArray::Scalar)1.3f);
        EXPECT_EQ(val2.y(), (typename TArray::Scalar)2.3f);

        // double setter and getter functions
        myArray->setValueDouble(0, as_double_t<typename TArray::Element>(11.1, 12.1));
        myArray->setValueDouble(1, as_double_t<typename TArray::Element>(11.2, 12.2));
        myArray->setValueDouble(2, as_double_t<typename TArray::Element>(11.3, 12.3));
        as_double_t<typename TArray::Element> val0d = myArray->getValueDouble(0);
        as_double_t<typename TArray::Element> val1d = myArray->getValueDouble(1);
        as_double_t<typename TArray::Element> val2d = myArray->getValueDouble(2);
        EXPECT_EQ(val0d.x(), (typename TArray::Scalar)11.1);
        EXPECT_EQ(val0d.y(), (typename TArray::Scalar)12.1);
        EXPECT_EQ(val1d.x(), (typename TArray::Scalar)11.2);
        EXPECT_EQ(val1d.y(), (typename TArray::Scalar)12.2);
        EXPECT_EQ(val2d.x(), (typename TArray::Scalar)11.3);
        EXPECT_EQ(val2d.y(), (typename TArray::Scalar)12.3);

        // int64_t setter and getter functions
        myArray->setValueInt64(0, as_int64_t<typename TArray::Element>(2, 3));
        myArray->setValueInt64(1, as_int64_t<typename TArray::Element>(3, 4));
        myArray->setValueInt64(2, as_int64_t<typename TArray::Element>(4, 5));
        as_int64_t<typename TArray::Element> val0i = myArray->getValueInt64(0);
        as_int64_t<typename TArray::Element> val1i = myArray->getValueInt64(1);
        as_int64_t<typename TArray::Element> val2i = myArray->getValueInt64(2);
        EXPECT_EQ(val0i.x(), (typename TArray::Scalar)2);
        EXPECT_EQ(val0i.y(), (typename TArray::Scalar)3);
        EXPECT_EQ(val1i.x(), (typename TArray::Scalar)3);
        EXPECT_EQ(val1i.y(), (typename TArray::Scalar)4);
        EXPECT_EQ(val2i.x(), (typename TArray::Scalar)4);
        EXPECT_EQ(val2i.y(), 5);

        // resize
        myArray->setSize(4);
        EXPECT_EQ(myArray->getSize(), 4);
    }

    /**
     * @brief Test function for arrays with three components.
     * @tparam TArray Array type to run the tests for.
     */
    template <typename TArray>
    void test_Array3()
    {
        auto myArray               = std::make_unique<TArray>();
        Eigen::Index numComponents = myArray->getNumberOfComponents();
        EXPECT_EQ(numComponents, 3);

        myArray->setSize(3);
        EXPECT_EQ(myArray->getSize(), 3);

        myArray->setName("test string");
        EXPECT_EQ(myArray->getName(), "test string");

        // regular setter and getter functions
        myArray->setValue(0, typename TArray::Element(1.1f, 2.1f, 3.1f));
        myArray->setValue(1, typename TArray::Element(1.2f, 2.2f, 3.2f));
        myArray->setValue(2, typename TArray::Element(1.3f, 2.3f, 3.3f));
        typename TArray::Element val0 = myArray->getValue(0);
        typename TArray::Element val1 = myArray->getValue(1);
        typename TArray::Element val2 = myArray->getValue(2);
        EXPECT_EQ(val0.x(), (typename TArray::Scalar)1.1f);
        EXPECT_EQ(val0.y(), (typename TArray::Scalar)2.1f);
        EXPECT_EQ(val0.z(), (typename TArray::Scalar)3.1f);
        EXPECT_EQ(val1.x(), (typename TArray::Scalar)1.2f);
        EXPECT_EQ(val1.y(), (typename TArray::Scalar)2.2f);
        EXPECT_EQ(val1.z(), (typename TArray::Scalar)3.2f);
        EXPECT_EQ(val2.x(), (typename TArray::Scalar)1.3f);
        EXPECT_EQ(val2.y(), (typename TArray::Scalar)2.3f);
        EXPECT_EQ(val2.z(), (typename TArray::Scalar)3.3f);

        // double setter and getter functions
        myArray->setValueDouble(0, as_double_t<typename TArray::Element>(11.1, 12.1, 13.1));
        myArray->setValueDouble(1, as_double_t<typename TArray::Element>(11.2, 12.2, 13.2));
        myArray->setValueDouble(2, as_double_t<typename TArray::Element>(11.3, 12.3, 13.3));
        as_double_t<typename TArray::Element> val0d = myArray->getValueDouble(0);
        as_double_t<typename TArray::Element> val1d = myArray->getValueDouble(1);
        as_double_t<typename TArray::Element> val2d = myArray->getValueDouble(2);
        EXPECT_EQ(val0d.x(), (typename TArray::Scalar)11.1);
        EXPECT_EQ(val0d.y(), (typename TArray::Scalar)12.1);
        EXPECT_EQ(val0d.z(), (typename TArray::Scalar)13.1);
        EXPECT_EQ(val1d.x(), (typename TArray::Scalar)11.2);
        EXPECT_EQ(val1d.y(), (typename TArray::Scalar)12.2);
        EXPECT_EQ(val1d.z(), (typename TArray::Scalar)13.2);
        EXPECT_EQ(val2d.x(), (typename TArray::Scalar)11.3);
        EXPECT_EQ(val2d.y(), (typename TArray::Scalar)12.3);
        EXPECT_EQ(val2d.z(), (typename TArray::Scalar)13.3);

        // int64_t setter and getter functions
        myArray->setValueInt64(0, as_int64_t<typename TArray::Element>(2, 3, 4));
        myArray->setValueInt64(1, as_int64_t<typename TArray::Element>(3, 4, 5));
        myArray->setValueInt64(2, as_int64_t<typename TArray::Element>(4, 5, 6));
        as_int64_t<typename TArray::Element> val0i = myArray->getValueInt64(0);
        as_int64_t<typename TArray::Element> val1i = myArray->getValueInt64(1);
        as_int64_t<typename TArray::Element> val2i = myArray->getValueInt64(2);
        EXPECT_EQ(val0i.x(), (typename TArray::Scalar)2);
        EXPECT_EQ(val0i.y(), (typename TArray::Scalar)3);
        EXPECT_EQ(val0i.z(), (typename TArray::Scalar)4);
        EXPECT_EQ(val1i.x(), (typename TArray::Scalar)3);
        EXPECT_EQ(val1i.y(), (typename TArray::Scalar)4);
        EXPECT_EQ(val1i.z(), (typename TArray::Scalar)5);
        EXPECT_EQ(val2i.x(), (typename TArray::Scalar)4);
        EXPECT_EQ(val2i.y(), (typename TArray::Scalar)5);
        EXPECT_EQ(val2i.z(), (typename TArray::Scalar)6);

        // resize
        myArray->setSize(4);
        EXPECT_EQ(myArray->getSize(), 4);
    }

    /**
     * @brief Test function for arrays with four components.
     * @tparam TArray Array type to run the tests for.
     */
    template <typename TArray>
    void test_Array4()
    {
        auto myArray               = std::make_unique<TArray>();
        Eigen::Index numComponents = myArray->getNumberOfComponents();
        EXPECT_EQ(numComponents, 4);

        myArray->setSize(3);
        EXPECT_EQ(myArray->getSize(), 3);

        myArray->setName("test string");
        EXPECT_EQ(myArray->getName(), "test string");

        // regular setter and getter functions
        myArray->setValue(0, typename TArray::Element(1.1f, 2.1f, 3.1f, 4.1f));
        myArray->setValue(1, typename TArray::Element(1.2f, 2.2f, 3.2f, 4.2f));
        myArray->setValue(2, typename TArray::Element(1.3f, 2.3f, 3.3f, 4.3f));
        typename TArray::Element val0 = myArray->getValue(0);
        typename TArray::Element val1 = myArray->getValue(1);
        typename TArray::Element val2 = myArray->getValue(2);
        EXPECT_EQ(val0.x(), (typename TArray::Scalar)1.1f);
        EXPECT_EQ(val0.y(), (typename TArray::Scalar)2.1f);
        EXPECT_EQ(val0.z(), (typename TArray::Scalar)3.1f);
        EXPECT_EQ(val0.w(), (typename TArray::Scalar)4.1f);
        EXPECT_EQ(val1.x(), (typename TArray::Scalar)1.2f);
        EXPECT_EQ(val1.y(), (typename TArray::Scalar)2.2f);
        EXPECT_EQ(val1.z(), (typename TArray::Scalar)3.2f);
        EXPECT_EQ(val1.w(), (typename TArray::Scalar)4.2f);
        EXPECT_EQ(val2.x(), (typename TArray::Scalar)1.3f);
        EXPECT_EQ(val2.y(), (typename TArray::Scalar)2.3f);
        EXPECT_EQ(val2.z(), (typename TArray::Scalar)3.3f);
        EXPECT_EQ(val2.w(), (typename TArray::Scalar)4.3f);

        // double setter and getter functions
        myArray->setValueDouble(0, as_double_t<typename TArray::Element>(11.1, 12.1, 13.1, 14.1));
        myArray->setValueDouble(1, as_double_t<typename TArray::Element>(11.2, 12.2, 13.2, 14.2));
        myArray->setValueDouble(2, as_double_t<typename TArray::Element>(11.3, 12.3, 13.3, 14.3));
        as_double_t<typename TArray::Element> val0d = myArray->getValueDouble(0);
        as_double_t<typename TArray::Element> val1d = myArray->getValueDouble(1);
        as_double_t<typename TArray::Element> val2d = myArray->getValueDouble(2);
        EXPECT_EQ(val0d.x(), (typename TArray::Scalar)11.1);
        EXPECT_EQ(val0d.y(), (typename TArray::Scalar)12.1);
        EXPECT_EQ(val0d.z(), (typename TArray::Scalar)13.1);
        EXPECT_EQ(val0d.w(), (typename TArray::Scalar)14.1);
        EXPECT_EQ(val1d.x(), (typename TArray::Scalar)11.2);
        EXPECT_EQ(val1d.y(), (typename TArray::Scalar)12.2);
        EXPECT_EQ(val1d.z(), (typename TArray::Scalar)13.2);
        EXPECT_EQ(val1d.w(), (typename TArray::Scalar)14.2);
        EXPECT_EQ(val2d.x(), (typename TArray::Scalar)11.3);
        EXPECT_EQ(val2d.y(), (typename TArray::Scalar)12.3);
        EXPECT_EQ(val2d.z(), (typename TArray::Scalar)13.3);
        EXPECT_EQ(val2d.w(), (typename TArray::Scalar)14.3);

        // int64_t setter and getter functions
        myArray->setValueInt64(0, as_int64_t<typename TArray::Element>(2, 3, 4, 5));
        myArray->setValueInt64(1, as_int64_t<typename TArray::Element>(3, 4, 5, 6));
        myArray->setValueInt64(2, as_int64_t<typename TArray::Element>(4, 5, 6, 7));
        as_int64_t<typename TArray::Element> val0i = myArray->getValueInt64(0);
        as_int64_t<typename TArray::Element> val1i = myArray->getValueInt64(1);
        as_int64_t<typename TArray::Element> val2i = myArray->getValueInt64(2);
        EXPECT_EQ(val0i.x(), (typename TArray::Scalar)2);
        EXPECT_EQ(val0i.y(), (typename TArray::Scalar)3);
        EXPECT_EQ(val0i.z(), (typename TArray::Scalar)4);
        EXPECT_EQ(val0i.w(), (typename TArray::Scalar)5);
        EXPECT_EQ(val1i.x(), (typename TArray::Scalar)3);
        EXPECT_EQ(val1i.y(), (typename TArray::Scalar)4);
        EXPECT_EQ(val1i.z(), (typename TArray::Scalar)5);
        EXPECT_EQ(val1i.w(), (typename TArray::Scalar)6);
        EXPECT_EQ(val2i.x(), (typename TArray::Scalar)4);
        EXPECT_EQ(val2i.y(), (typename TArray::Scalar)5);
        EXPECT_EQ(val2i.z(), (typename TArray::Scalar)6);
        EXPECT_EQ(val2i.w(), (typename TArray::Scalar)7);

        // resize
        myArray->setSize(4);
        EXPECT_EQ(myArray->getSize(), 4);
    }

    /**
     * @brief Tests that perform basic operations on an array.
     */
    void test_Array_manipulation()
    {
        // IsEqual
        auto arrayA = std::make_unique<Array1i>();
        arrayA->setValues({ 1, 3, 2 });
        auto arrayA_cpy = std::make_unique<Array1i>();
        arrayA_cpy->setValues({ 1, 3, 2 });
        EXPECT_TRUE(arrayA->isEqual(arrayA_cpy.get()));

        // Append
        auto arrayB = std::make_unique<Array1i>();
        arrayB->setValues({ 5, 2, 4, 6, 1 });
        arrayA->append(arrayB.get());
        auto array_append_gt = std::make_unique<Array1i>();
        array_append_gt->setValues({ 1, 3, 2, 5, 2, 4, 6, 1 });
        EXPECT_TRUE(arrayA->isEqual(array_append_gt.get()));

        // Prepend
        auto arrayB2 = std::make_unique<Array1i>();
        arrayB2->setValues({ 3, 2, 3, 6 });
        arrayA->prepend(arrayB2.get());
        auto array_prepend_gt = std::make_unique<Array1i>();
        array_prepend_gt->setValues({ 3, 2, 3, 6, 1, 3, 2, 5, 2, 4, 6, 1 });
        EXPECT_TRUE(arrayA->isEqual(array_prepend_gt.get()));

        // Reverse
        arrayA->reverse();
        auto array_reverse_gt = std::make_unique<Array1i>();
        array_reverse_gt->setValues({ 1, 6, 4, 2, 5, 2, 3, 1, 6, 3, 2, 3 });
        EXPECT_TRUE(arrayA->isEqual(array_reverse_gt.get()));

        // Append vector
        std::shared_ptr<Array2d> arrayC = std::make_unique<Array2d>();
        arrayC->setValues({ Eigen::Vector2d(1, 1.1), Eigen::Vector2d(2, 2.2), Eigen::Vector2d(3, 3.3) });
        auto arrayD = std::make_unique<Array2d>();
        arrayD->setValues({ Eigen::Vector2d(4, 4.4), Eigen::Vector2d(5, 5.5), Eigen::Vector2d(6, 6.6), Eigen::Vector2d(7, 7.7) });
        auto arrayE = arrayC;
        arrayE->append(arrayD.get());
        auto arrayE_gt = std::make_unique<Array2d>();
        arrayE_gt->setValues({ Eigen::Vector2d(1, 1.1), Eigen::Vector2d(2, 2.2), Eigen::Vector2d(3, 3.3), Eigen::Vector2d(4, 4.4), Eigen::Vector2d(5, 5.5), Eigen::Vector2d(6, 6.6), Eigen::Vector2d(7, 7.7) });
        EXPECT_TRUE(arrayE->isEqual(arrayE_gt.get()));

        // First, Last
        EXPECT_EQ(arrayD->first(), Eigen::Vector2d(4, 4.4));
        EXPECT_EQ(arrayD->last(), Eigen::Vector2d(7, 7.7));

        // GetNumberOfComponents
        EXPECT_EQ(arrayA->getNumberOfComponents(), 1);
        EXPECT_EQ(arrayD->getNumberOfComponents(), 2);

        // Min, Max
        auto arrayF = std::make_unique<Array1i>();
        arrayF->setValues({ 3, 1, 6, 4, 2, 7 });
        EXPECT_EQ(arrayF->getMin().x(), 1);
        EXPECT_EQ(arrayF->getMax().x(), 7);

        // RemoveLast
        auto arrayG = std::make_unique<Array2d>();
        arrayG->setValues({ Eigen::Vector2d(1, 1.1), Eigen::Vector2d(2, 2.2), Eigen::Vector2d(3, 3.3), Eigen::Vector2d(4, 4.4), Eigen::Vector2d(5, 5.5), Eigen::Vector2d(6, 6.6), Eigen::Vector2d(7, 7.7) });
        arrayG->removeLast(5);
        auto arrayG_gt = std::make_unique<Array2d>();
        arrayG_gt->setValues({ Eigen::Vector2d(1, 1.1), Eigen::Vector2d(2, 2.2) });
        EXPECT_TRUE(arrayG->isEqual(arrayG_gt.get()));

        // RemoveFirst
        auto arrayH = std::make_unique<Array2d>();
        arrayH->setValues({ Eigen::Vector2d(1, 1.1), Eigen::Vector2d(2, 2.2), Eigen::Vector2d(3, 3.3), Eigen::Vector2d(4, 4.4), Eigen::Vector2d(5, 5.5), Eigen::Vector2d(6, 6.6), Eigen::Vector2d(7, 7.7) });
        arrayH->removeFirst(5);
        auto arrayH_gt = std::make_unique<Array2d>();
        arrayH_gt->setValues({ Eigen::Vector2d(6, 6.6), Eigen::Vector2d(7, 7.7) });
        EXPECT_TRUE(arrayH->isEqual(arrayH_gt.get()));

        // GetElementSize, GetSize, GetSizeInBytes
        EXPECT_EQ(arrayC->getElementSizeInBytes(), 16);
        EXPECT_EQ(arrayC->getSize(), 7);
        EXPECT_EQ(arrayC->getSizeInBytes(), 112);

        // Resize
        auto arrayI = std::make_unique<Array2d>();
        arrayI->setValues({ Eigen::Vector2d(1, 1.1), Eigen::Vector2d(2, 2.2), Eigen::Vector2d(3, 3.3) });
        arrayI->setSize(4);
        arrayI->setValue(3, Eigen::Vector2d(0, 0));
        auto arrayI_gt = std::make_unique<Array2d>();
        arrayI_gt->setValues({ Eigen::Vector2d(1, 1.1), Eigen::Vector2d(2, 2.2), Eigen::Vector2d(3, 3.3), Eigen::Vector2d(0, 0) });
        EXPECT_TRUE(arrayI->isEqual(arrayI_gt.get()));
        arrayI->setSize(2);
        auto arrayI_gt2 = std::make_unique<Array2d>();
        arrayI_gt2->setValues({ Eigen::Vector2d(1, 1.1), Eigen::Vector2d(2, 2.2) });
        EXPECT_TRUE(arrayI->isEqual(arrayI_gt2.get()));

        // SortAscending, SortDescending
        auto arrayJ = std::make_unique<Array1i>();
        arrayJ->setValues({ 3, 1, 6, 4, 2, 7 });
        arrayJ->sortAscending();
        auto arrayJ_ascend = std::make_unique<Array1i>();
        arrayJ_ascend->setValues({ 1, 2, 3, 4, 6, 7 });
        EXPECT_TRUE(arrayJ->isEqual(arrayJ_ascend.get()));
        arrayJ->sortDescending();
        auto arrayJ_descend = std::make_unique<Array1i>();
        arrayJ_descend->setValues({ 7, 6, 4, 3, 2, 1 });
        EXPECT_TRUE(arrayJ->isEqual(arrayJ_descend.get()));

        // Clone
        auto arrayK = std::make_unique<Array1i>();
        arrayK->setValues({ 3, 1, 6, 4, 2, 7 });
        auto cloned       = std::shared_ptr<Array1i>(arrayK->clone());
        auto arrayK_clone = std::dynamic_pointer_cast<Array1i, Object>(cloned);
        EXPECT_TRUE(arrayK->isEqual(arrayK_clone.get()));
        arrayK->setValue(0, 4);
        EXPECT_FALSE(arrayK->isEqual(arrayK_clone.get()));
        arrayK_clone->setValue(0, 4);
        EXPECT_TRUE(arrayK->isEqual(arrayK_clone.get()));
    }

    TEST(core, array)
    {
        test_ArrayX<ArrayXf>();
        test_Array1<Array1f>();
        test_Array2<Array2f>();
        test_Array3<Array3f>();
        test_Array4<Array4f>();

        test_ArrayX<ArrayXd>();
        test_Array1<Array1d>();
        test_Array2<Array2d>();
        test_Array3<Array3d>();
        test_Array4<Array4d>();

        test_ArrayX<ArrayXi>();
        test_Array1<Array1i>();
        test_Array2<Array2i>();
        test_Array3<Array3i>();
        test_Array4<Array4i>();

        test_ArrayX<SparseArrayXf>();
        test_Array1<SparseArray1f>();
        test_Array2<SparseArray2f>();
        test_Array3<SparseArray3f>();
        test_Array4<SparseArray4f>();

        test_ArrayX<SparseArrayXd>();
        test_Array1<SparseArray1d>();
        test_Array2<SparseArray2d>();
        test_Array3<SparseArray3d>();
        test_Array4<SparseArray4d>();

        test_ArrayX<SparseArrayXi>();
        test_Array1<SparseArray1i>();
        test_Array2<SparseArray2i>();
        test_Array3<SparseArray3i>();
        test_Array4<SparseArray4i>();

        test_Array_manipulation();
    }
}
