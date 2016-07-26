#pragma once

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <functional>
#include <stdexcept>
#include <bitset>
#include <iostream>
#include <sstream>

#include "../value.h"
#include "../dhtrunner.h"
#include "../infohash.h"

namespace dht {
namespace indexation {

/*!
 * @class   Prefix
 * @brief   A blob structure which prefixes a Key in the PHT.
 * @details
 * Since the PHT structure is a "trie", every node in this structure have a
 * label which is defined by the path from the root of the trie to the node. If
 * the node in question is a leaf, *the label is a prefix of all the keys
 * contained in the leaf*.
 */
struct Prefix {
    Prefix() {}
    Prefix(InfoHash h) : size_(h.size() * 8), content_(h.begin(), h.end()) { }
    Prefix(const Blob& d, const Blob& f={}) : size_(d.size()*8), flags_(f), content_(d) { }

    Prefix(const Prefix& p, size_t first) :
        size_(std::min(first, p.content_.size()*8)),
        content_(Blob(p.content_.begin(), p.content_.begin()+size_/8))
    {
        auto rem = size_ % 8;
        if ( not flags_.empty() ) {
            flags_ = Blob(p.flags_.begin(), p.flags_.begin()+size_/8);
            if (rem)
                flags_.push_back(p.flags_[size_/8] & (0xFF << (8 - rem))); 
        }

        if (rem)
            content_.push_back(p.content_[size_/8] & (0xFF << (8 - rem)));
    }

    Prefix getPrefix(ssize_t len) const {
        if ((size_t)std::abs(len) > size_)
            throw std::out_of_range("len larger than prefix size.");
        if (len < 0)
            len += size_;
        return Prefix(*this, len);
    }

    /**
     * Method for getting the state of the bit at the position pos.
     * @param pos : Pos of the needed bit
     * @return : true if the bit is at 1
     *           false otherwise
     * @throw out_of_range Throw out of range if the bit at 'pos' does not exist
     */
    bool isFlagActive(size_t pos) const {
        return flags_.empty() or isActiveBit(flags_, pos);
    }

    bool isContentBitActive(size_t pos) const {
        return isActiveBit(content_, pos);
    }

    Prefix getFullSize() { return Prefix(*this, content_.size()*8); }

    /**
     * This methods gets the prefix of its sibling in the PHT structure.
     *
     * @return The prefix of this sibling.
     */
    Prefix getSibling() const {
        Prefix copy = *this;
        copy.swapContentBit(size_ - 1);

        return copy;
    }

    InfoHash hash() const {
        Blob copy(content_);
        copy.push_back(size_);
        return InfoHash::get(copy);
    }

    std::string toString() const {
        std::stringstream ss;

        ss << "Prefix : " << std::endl << "\tContent_ : ";
        ss << blobToString(content_);
        ss << std::endl;

        ss << "\tFlags_ :   ";
        ss << blobToString(flags_);
        ss << std::endl;

        return ss.str();
    }

    static inline unsigned commonBits(const Prefix& p1, const Prefix& p2) {
        unsigned i, j;
        uint8_t x;
        auto longest_prefix_size = std::min(p1.size_, p2.size_);

        for (i = 0; i < longest_prefix_size; i++) {
            if (p1.content_.data()[i] != p2.content_.data()[i]
                or not p1.isFlagActive(i) 
                or not p2.isFlagActive(i) ) {

                    break;
            }
        }

        if (i == longest_prefix_size)
            return 8*longest_prefix_size;

        x = p1.content_.data()[i] ^ p2.content_.data()[i];

        j = 0;
        while ((x & 0x80) == 0) {
            x <<= 1;
            j++;
        }

        return 8 * i + j;
    }

    /**
     * This method swap the bit a the position 'bit'
     *
     * @param bit Position of the bit to swap
     * @return The prefix with the bit at position 'bit' swapped
     * @throw out_of_range Throw out of range if bit does not exist
     */
    void swapContentBit(size_t bit) {
        swapBit(content_, bit);
    }

    void swapFlagBit(size_t bit) {
        swapBit(flags_, bit);
    }

    void addPaddingContent(size_t size) {
        content_ = addPadding(content_, size);
    }

    void updateFlags() {
        // Fill first known bit
        auto csize = size_ - flags_.size() * 8;
        while(csize >= 8) {
            flags_.push_back(0xFF);
            csize -= 8;
        }

        // if needed fill remaining bit
        if ( csize )
            flags_.push_back(0xFF << (8 - csize));

        // Complet vector space missing
        for ( auto i = flags_.size(); i < content_.size(); i++ )
            flags_.push_back(0xFF);
    }

    size_t size_ {0};

    Blob flags_ {};
    Blob content_ {};
private:
    std::string blobToString(const Blob &bl) const {
        std::stringstream ss;

        auto bn = size_ % 8;
        auto n = size_ / 8;

        for (size_t i = 0; i < bl.size(); i++)
            ss << std::bitset<8>(bl[i]) << " ";
        if (bn)
            for (unsigned b=0; b < bn; b++)
                ss << (char)((bl[n] & (1 << (7 - b))) ? '1':'0');

        return ss.str();
    }

    Blob addPadding(Blob toP, size_t size) {
        Blob copy = toP;
        for ( auto i = copy.size(); i < size; i++ )
            copy.push_back(0);

        swapBit(copy, size_ + 1);
        return copy;
    }

    bool isActiveBit(const Blob &b, size_t pos) const {
        if ( pos >= size_ )
            throw std::out_of_range("Can't detect active bit at pos, pos larger than prefix size or empty prefix");

        return ((b[pos / 8] >> (7 - (pos % 8)) ) & 1) == 1;
    }

    void swapBit(Blob &b, size_t bit) {
        if ( bit >= b.size() * 8 )
            throw std::out_of_range("bit larger than prefix size.");

        size_t offset_bit = (8 - bit) % 8;
        b[bit / 8] ^= (1 << offset_bit);
    }
};

using Value = std::pair<InfoHash, dht::Value::Id>;
struct IndexEntry : public dht::Value::Serializable<IndexEntry> {
    static const ValueType TYPE;

    virtual void unpackValue(const dht::Value& v) {
        Serializable<IndexEntry>::unpackValue(v);
        name = v.user_type;
    }

    virtual dht::Value packValue() const {
        auto pack = Serializable<IndexEntry>::packValue();
        pack.user_type = name;
        return pack;
    }

    Blob prefix;
    Value value;
    std::string name;
    MSGPACK_DEFINE_MAP(prefix, value);
};

class Pht {
    static constexpr const char* INVALID_KEY = "Key does not match the PHT key spec.";

    /* Prefixes the user_type for all dht values put on the DHT */
    static constexpr const char* INDEX_PREFIX = "index.pht.";

public:

    /* This is the maximum number of entries per node. This parameter is
     * critical and influences the traffic a lot during a lookup operation.
     */
    static constexpr const size_t MAX_NODE_ENTRY_COUNT {16};

    /* A key for a an index entry */
    using Key = std::map<std::string, Blob>;
    /* Specifications of the keys. It defines the number, the length and the
     * serialization order of fields. */
    using KeySpec = std::map<std::string, size_t>;

    using RealInsertCallback = std::function<void(std::shared_ptr<Prefix> p, IndexEntry entry )>;
    using LookupCallback = std::function<void(std::vector<std::shared_ptr<Value>>& values, Prefix p)>;
    using LookupCallbackWrapper = std::function<void(std::vector<std::shared_ptr<IndexEntry>>& values, Prefix p)>;

    typedef void (*LookupCallbackRaw)(std::vector<std::shared_ptr<Value>>* values, Prefix* p, void *user_data);
    static LookupCallback
    bindLookupCb(LookupCallbackRaw raw_cb, void* user_data) {
        if (not raw_cb) return {};
        return [=](std::vector<std::shared_ptr<Value>>& values, Prefix p) {
            raw_cb((std::vector<std::shared_ptr<Value>>*) &values, (Prefix*) &p, user_data);
        };
    }
    using LookupCallbackSimple = std::function<void(std::vector<std::shared_ptr<Value>>& values)>;
    typedef void (*LookupCallbackSimpleRaw)(std::vector<std::shared_ptr<Value>>* values, void *user_data);
    static LookupCallbackSimple
    bindLookupCbSimple(LookupCallbackSimpleRaw raw_cb, void* user_data) {
        if (not raw_cb) return {};
        return [=](std::vector<std::shared_ptr<Value>>& values) {
            raw_cb((std::vector<std::shared_ptr<Value>>*) &values, user_data);
        };
    }

    Pht(std::string name, KeySpec k_spec, std::shared_ptr<DhtRunner> dht)
        : name_(INDEX_PREFIX + name), canary_(name_ + ".canary"), keySpec_(k_spec), dht_(dht) {}

    virtual ~Pht () { }

    /**
     * Lookup a key for a value.
     */
    void lookup(Key k, LookupCallback cb = {}, DoneCallbackSimple done_cb = {}, bool exact_match = true);
    void lookup(Key k, LookupCallbackSimple cb = {}, DoneCallbackSimple done_cb = {}, bool exact_match = true)
    {
        lookup(k, [=](std::vector<std::shared_ptr<Value>>& values, Prefix) { cb(values); }, done_cb, exact_match);
    }

    /**
     * Adds an entry into the index.
     */
    void insert(Key k, Value v, DoneCallbackSimple done_cb = {}) {
        Prefix p = linearize(k);

        auto lo = std::make_shared<int>(0);
        auto hi = std::make_shared<int>(p.size_);

        IndexEntry entry;
        entry.value = v;
        entry.prefix = p.content_;
        entry.name = name_;

        Pht::insert(p, entry, lo, hi, clock::now(), done_cb);
    }

private:
    void insert( Prefix kp, IndexEntry entry, std::shared_ptr<int> lo, std::shared_ptr<int> hi, time_point time_p,
                 DoneCallbackSimple done_cb = {});

    class Cache {
    public:
        /**
         * Insert all needed node into the tree according to a prefix
         * @param p : Prefix that we need to insert
         */
        void insert(const Prefix& p);

        /**
         * Lookup into the tree to return the maximum prefix length in the cache tree
         *
         * @param p : Prefix that we are looking for
         * @return  : The size of the longest prefix known in the cache between 0 and p.size_
         */
        int lookup(const Prefix& p);

    private:
        static constexpr const size_t MAX_ELEMENT {1024};
        static constexpr const std::chrono::minutes NODE_EXPIRE_TIME {5};

        struct Node {
            time_point last_reply;           /* Made the assocation between leaves and leaves multimap */
            std::shared_ptr<Node> parent;    /* Share_ptr to the parent, it allow the self destruction of tree */
            std::weak_ptr<Node> left_child;  /* Left child, for bit equal to 1 */
            std::weak_ptr<Node> right_child; /* Right child, for bit equal to 0 */
        };

        std::weak_ptr<Node> root_;                         /* Root of the tree */

        /**
         * This mutlimap contains all prefix insert in the tree in time order
         * We could then delete the last one if there is too much node
         * The tree will self destroy is branch ( thanks to share_ptr )
         */
        std::multimap<time_point, std::shared_ptr<Node>> leaves_;
    };

    /**
     * Performs a step in the lookup operation. Each steps are performed
     * asynchronously.
     */
    void lookupStep(Prefix k, std::shared_ptr<int> lo, std::shared_ptr<int> hi,
            std::shared_ptr<std::vector<std::shared_ptr<IndexEntry>>> vals, LookupCallbackWrapper cb,
            DoneCallbackSimple done_cb, std::shared_ptr<unsigned> max_common_prefix_len,
            int start = -1, bool all_values = false);
    
    Prefix zcurve(const std::vector<Prefix>& all_prefix) const {
        Prefix p;
        if ( all_prefix.size() == 1 ) return all_prefix[0];

        for ( size_t j = 0, bit = 0; j < all_prefix[0].content_.size(); j++) {
            uint8_t mask = 0x80;
            for ( int i = 0; i < 8; ) {
                uint8_t content = 0;
                uint8_t flags = 0;
                for ( int k = 0 ; k < 8; k++, bit++ ) {
                    auto diff = k - i;
                    auto x = all_prefix[bit].content_[j] & mask;
                    auto y = all_prefix[bit].flags_[j] & mask;
                    content |= ( diff >= 0 ) ? x >> diff : x << std::abs(diff);
                    flags   |= ( diff >= 0 ) ? y >> diff : y << std::abs(diff);

                    if ( bit == all_prefix.size() - 1 ) { bit = -1; ++i; mask >>= 1; }
                }
                p.content_.push_back(content);
                p.flags_.push_back(flags);
                p.size_ += 8;
            }
        }
        std::cerr << p.toString() << std::endl;
        return p;   
    }

    /**
     * Linearizes the key into a unidimensional key. A pht only takes
     * unidimensional key.
     *
     * @param Key  The initial key.
     *
     * @return the prefix of the linearized key.
     */
    virtual Prefix linearize(Key k) const {
        if (not validKey(k)) { throw std::invalid_argument(INVALID_KEY); }

        std::vector<Prefix> all_prefix;
        auto max = std::max_element(keySpec_.begin(), keySpec_.end(), 
            [&](const std::pair<std::string, size_t>& a, const std::pair<std::string, size_t>& b) {
                return a.second < b.second;
            })->second + 1;

        for ( auto const& it : k ) {
            Prefix p = Blob {it.second.begin(), it.second.end()};
            p.addPaddingContent(max);
            p.updateFlags();
            all_prefix.push_back(p);

            std::cerr << p.toString() << std::endl;
        }

        return zcurve(all_prefix);
    };

    /**
     * Tells if the key is valid according to the key spec.
     */
    bool validKey(const Key& k) const {
        return k.size() == keySpec_.size() and
            std::equal(k.begin(), k.end(), keySpec_.begin(),
                [&](const Key::value_type& key, const KeySpec::value_type& key_spec) {
                    return key.first == key_spec.first and key.second.size() <= key_spec.second;
                }
            );
    }

    /**
     * Looking where to put the data cause if there i free space on the node above then this node will became the real leave.
     *
     * @param p       Share_ptr on the Prefix to check
     * @param entry   The entry to put at the prefix p
     * @param end_cb  Callback to use at the end of counting
     */
    void getRealPrefix(std::shared_ptr<Prefix> p, IndexEntry entry, RealInsertCallback end_cb ) {

        if ( p->size_ == 0 ) {
            end_cb(p, std::move(entry));
            return;
        }

        auto total = std::make_shared<unsigned int>(0); /* Will contains the total number of data on 3 nodes */
        auto ended = std::make_shared<unsigned int>(0); /* Just indicate how many have end */

        auto parent = std::make_shared<Prefix>(p->getPrefix(-1));
        auto sibling = std::make_shared<Prefix>(p->getSibling());

        auto pht_filter = [&](const dht::Value& v) {
            return v.user_type.compare(0, name_.size(), name_) == 0;
        };

        /* Lambda will count total number of data node */
        auto count = [=]( const std::shared_ptr<dht::Value> value ) {
            if ( value->user_type != canary_)
                (*total)++;

            return true;
        };

        auto on_done = [=] ( bool ) {
            (*ended)++;
            /* Only the last one do the CallBack*/
            if  ( *ended == 3 ) {
                if ( *total < MAX_NODE_ENTRY_COUNT )
                    end_cb(parent, std::move(entry));
                else
                    end_cb(p, std::move(entry));
            }
        };

        dht_->get(parent->hash(),
            count,
            on_done,
            pht_filter
        );

        dht_->get(p->hash(),
            count,
            on_done,
            pht_filter
        );

        dht_->get(sibling->hash(),
            count,
            on_done,
            pht_filter
        );
    }

    /**
     * Tells if the key is valid according to the key spec.
     */
    void checkPhtUpdate(Prefix p, IndexEntry entry, time_point time_p) {

        Prefix full = entry.prefix;
        if ( p.size_ >= full.size_ ) return;

        auto next_prefix = full.getPrefix( p.size_ + 1 ); 

        dht_->listen(next_prefix.hash(),
            [=](const std::shared_ptr<dht::Value> &value) {
                if (value->user_type == canary_) {
                    insert(p, entry, std::make_shared<int>(0), std::make_shared<int>(p.size_), time_p, nullptr);

                    /* Cancel listen since we found where we need to update*/
                    return false;
                }

                return true;
            },
            [=](const dht::Value& v) {
                /* Filter value v thats start with the same name as ours */
                return v.user_type.compare(0, name_.size(), name_) == 0;
            }
        );
    }

    size_t foundSplitLocation(Prefix compared, std::shared_ptr<std::vector<std::shared_ptr<IndexEntry>>> vals) {
        for ( size_t i = 0; i < compared.size_; i++ ) 
            for ( auto const& v : *vals)
                if ( Prefix(v->prefix).isContentBitActive(i) != compared.isContentBitActive(i) )
                    return i;

        return compared.size_ - 1;
    }

    void split(Prefix insert, std::shared_ptr<std::vector<std::shared_ptr<IndexEntry>>> vals, IndexEntry entry, RealInsertCallback end_cb ) {
        auto full = Prefix(entry.prefix);

        auto loc = foundSplitLocation(full, vals);
        auto prefix_to_insert = std::make_shared<Prefix>(full.getPrefix(loc + 1));

        for (; loc > insert.size_; --loc)
            updateCanary(full.getPrefix(loc));
        
        end_cb(prefix_to_insert, entry);
}
    /**
     * Updates the canary token on the node responsible for the specified
     * Prefix.
     */
    void updateCanary(Prefix p);

    const std::string name_;
    const std::string canary_;
    const KeySpec keySpec_;
    Cache cache_;
    std::shared_ptr<DhtRunner> dht_;
};

} /* indexation  */
} /* dht */

