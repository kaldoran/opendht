#include <iostream>
#include <vector>
#include <cmath>
#include <bitset>
#include <string>
#include <sstream>

using Blob = std::vector<uint8_t>;

struct Prefix;

struct InfoPrefix {
    InfoPrefix() {}
    InfoPrefix(const Blob& d) : size_(d.size()*8), content_(d) {}
    InfoPrefix(const InfoPrefix& p, size_t first) :
        size_(std::min(first, p.content_.size()*8)),
        content_(Blob(p.content_.begin(), p.content_.begin()+size_/8))
    {
        auto rem = size_ % 8;
        if (rem)
            content_.push_back(p.content_[size_/8] & (0xFF << (8 - rem)));
    }

	Prefix toPrefix();

    InfoPrefix getPrefix(ssize_t len) const {
        if ((size_t)std::abs(len) > size_)
            throw std::out_of_range("len larger than prefix size.");
        if (len < 0)
            len += size_;
        return InfoPrefix(*this, len);
    }

    /**
     * This method swap the bit a the position 'bit' and return the new prefix
     * @param bit Position of the bit to swap
     * @return The prefix with the bit at position 'bit' swapped
     * @throw out_of_range Throw out of range if bit does not exist
     */
    virtual InfoPrefix swapBit(size_t bit) const {
        if ( bit >= size_ )
            throw std::out_of_range("bit larger than prefix size.");

        InfoPrefix copy = *this;

        size_t offset_bit = (8 - bit) % 8;
        copy.content_[bit / 8] ^= (1 << offset_bit);

        return copy;
    }

    /**
     * Method for getting the state of the bit at the position pos.
     * @param pos : Pos of the needed bit
     * @return : true if the bit is at 1
     *           false otherwise
     * @throw out_of_range Throw out of range if the bit at 'pos' does not exist
     */
    virtual bool isActiveBit(size_t pos) const {
        if ( pos >= size_ )
            throw std::out_of_range("Can't detect active bit at pos, pos larger than prefix size or empty prefix");

        return ((this->content_[pos / 8] >> (7 - (pos % 8)) ) & 1) == 1;
    }

    std::string toString() const {
        std::stringstream ss;
        auto bn = size_ % 8;
        auto n = size_ / 8;
        for (size_t i = 0; i<n; i++)
            ss << std::bitset<8>(content_[i]);
        if (bn)
            for (unsigned b=0; b<bn; b++)
                ss << (char)((content_[n] & (1 << (7 - b))) ? '1':'0');
        return ss.str();
    }

    InfoPrefix getFullSize() { return InfoPrefix(*this, content_.size()*8); }

	// Unexpect behavior if a prefix without padding is given.
	// No way to check if there is a prefix or if there is not
	InfoPrefix realPrefix() {
		for ( size_t i = content_.size() - 1; i >= 0; i-- ) { 
			uint8_t j = 0x01;
			for ( size_t pos = 0; pos < 8; pos++, j <<= 1 )
				if ( (content_[i] & j) != 0x00 )
					return getPrefix(pos + ( content_.size() - 1 - i) * 8);
		}
	 
		// Should never append
		return *this;
	}

    size_t size_ {0};
    Blob content_ {};
};

/*!
 * @class   Prefix
 * @brief   A blob structure which prefixes a Key in the PHT.
 * @details
 * Since the PHT structure is a "trie", every node in this structure have a
 * label which is defined by the path from the root of the trie to the node. If
 * the node in question is a leaf, *the label is a prefix of all the keys
 * contained in the leaf*.
 */
struct Prefix : InfoPrefix {
    Prefix() : InfoPrefix() {}
    // Prefix(InfoHash h) : size_(h.size() * 8), content_(h.begin(), h.end()) {}
    Prefix(const Blob& d) : InfoPrefix(d) {}
    Prefix(InfoPrefix const& infop) : InfoPrefix(infop) {}
    Prefix(const Prefix& p, size_t first) : InfoPrefix(p, first) {}

    /**
     * This methods gets the prefix of its sibling in the PHT structure.
     *
     * @return The prefix of this sibling.
     */
    Prefix getSibling() const {
        return swapBit(size_ - 1);
    }
/*
    InfoHash hash() const {
        Blob copy(content_);
        copy.push_back(size_);
        return InfoHash::get(copy);
    }
*/
    static inline unsigned commonBits(const Prefix& p1, const Prefix& p2) {
        unsigned i, j;
        uint8_t x;
        auto longest_prefix_size = std::min(p1.size_, p2.size_);

        for (i = 0; i < longest_prefix_size; i++) {
            if (p1.content_.data()[i] != p2.content_.data()[i])
                break;
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

    InfoPrefix toInfoPrefix() {

		Blob out(2 * content_.size()); // ceil(p.size_ * 2. / 8.);

		for ( size_t i = 0; i < content_.size(); i++ ) {
	    	for ( int b = 0, shift = 0; b < 8; b++, shift += 2) {
				if ( shift >= 8 ) 
					shift = 0;
	 
				/* 1 -> 11 [0xC0], 0 -> 10 [0x80], shift in right place */
	    		uint8_t mout = (( content_[i] & (0x80 >> b ) ) ? 0xC0 : 0x80) >> shift;
				out[i * 2 + (b >= 4)] |= mout; 
	    	}
		}

		return InfoPrefix(out);
	}

};


Prefix InfoPrefix::toPrefix() {
	Blob out(ceil(content_.size() / 2.)); // Output is on 1 bit instead of 2 for each data

	for ( size_t i = 0; i < content_.size(); i++ ) {
    	if ( (content_[i] | 0x55) != 0xFF )
    		throw std::domain_error("Can not convert to prefix with undefine bits"); 
 
    	uint8_t mask = 0x40;           // Mask for extracting bit a pos x               
    	int shift = ( i % 2 ) ? 3 : 1; // For first uin8_t need to move shift bit to left
    	                               // and to right for second uint8_t
    							       // First : 0101 0000 - Second : 1111 0111
    	                               // First : Bit pos 1 move to pos 0 - Second : Bit pos 1 move to pos 4
 
    	for ( size_t b = 1; b < 8; b += 2, mask >>= 2) { // for each bit of the uint8_t ( jump 1 of 2) and shift j mask
    		out[i / 2] ^= (
    						( i % 2 ) ? 
    							( content_[i] & mask ) >> shift  // Get bit pos l, move it to position by shifting of k
    					  	  : ( content_[i] & mask ) << shift
    					  );
 
    		( i % 2 ) ? --shift : ++shift; // Next shift to do
    	}
    }
 
    return out;
}


// http://stackoverflow.com/questions/29547628/convert-uint8-t-hex-value-to-binary
void toBinary(uint8_t a) {
	int tmp = 0;
    for(uint8_t i=0x800x80;i!=0;i>>=1) {
    	if ( tmp++ % 4 == 0 ) 
    		std::cout << " ";

        std::cout << ((a&i)?'1':'0');
    }
}

void toBinary(Blob p) {
	for (size_t i = 0; i < p.size(); i++ )
		toBinary(p[i]);
}

int main() {
	
	Blob p; // Test Blob	
	p.push_back(0xBA); // 1011 1010
	p.push_back(0xEA); // 1110 1010

	InfoPrefix pr(p);
	std::cout << "InfoPrefix normal : " << pr.toString() << std::endl;

	Prefix pp = pr.toPrefix();
	std::cout << "InfoPrefixrefix to Prefix : " << pp.toString() << std::endl;
	std::cout << "Prefix to InfoPrefix : " << pp.toInfoPrefix().toString() << std::endl;

    std::cout << std::endl;
    
    exit(0);
}


