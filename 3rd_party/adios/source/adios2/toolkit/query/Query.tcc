namespace adios2
{
namespace query
{

template <class T>
bool Range::CheckInterval(T &min, T &max) const
{
    bool isHit = false;
    std::stringstream convert(m_StrValue);
    T value;
    convert >> value;

    switch (m_Op)
    {
    case adios2::query::Op::GT:
        isHit = (max > value);
        break;
    case adios2::query::Op::LT:
        isHit = (min < value);
        break;
    case adios2::query::Op::GE:
        isHit = (max >= value);
        break;
    case adios2::query::Op::LE:
        isHit = (min <= value);
        break;
    case adios2::query::Op::EQ:
        isHit = (max >= value) && (min <= value);
        break;
    case adios2::query::Op::NE:
        isHit = !((max == value) && (min == value));
        break;
    default:
        break;
    }
    return isHit;
}

template <class T>
bool RangeTree::CheckInterval(T &min, T &max) const
{
    if (adios2::query::Relation::AND == m_Relation)
    {
        for (auto &range : m_Leaves)
            if (!range.CheckInterval(min, max))
                return false;

        for (auto &node : m_SubNodes)
            if (!node.CheckInterval(min, max))
                return false;

        return true; // even if no leaves or nodes
    }

    if (adios2::query::Relation::OR == m_Relation)
    {
        for (auto &range : m_Leaves)
            if (range.CheckInterval(min, max))
                return true;

        for (auto &node : m_SubNodes)
            if (node.CheckInterval(min, max))
                return true;

        return false;
    }

    // anything else are false
    return false;
}
}
}
