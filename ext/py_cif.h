#pragma once

#include "./py_utils.h"
#include <cctype>

namespace jnc {
namespace bio {

class CifParser {
public:
  std::deque<std::string> m_attributes;
  std::string m_category;

  std::istream &io;

  enum fsm_status {
    ST_START,
    ST_START_,
    ST_COMMENT,
    ST_ES,
    ST_SS,
    ST_SS_START,
    ST_SS_ESCAPE,
    ST_DS,
    ST_DS_START,
    ST_DS_ESCAPE,
    ST_MS,
    ST_MS_NL,
    ST_MS_START,
    ST_MS_ESCAPE,
    ST_CATEGORY_START,
    ST_CATEGORY,
    ST_ATTRIBUTE_START,
    ST_ATTRIBUTE,
    ST_VALUE,
    ST_ERROR,
    ST_EOF
  };

  std::map<fsm_status, std::function<fsm_status(char)>> fsm;

  enum word_t {
    WD_CATEGORY,
    WD_ATTRIBUTE,
    WD_VALUE,
  };

  // word parsing results
  struct word_parsing_t {
    std::string word;
    fsm_status status = ST_START;
    fsm_status pre_status = ST_START;
    word_t word_type;
    bool empty = true;
    int iline = 0;
    int col = 0;
  };

  word_parsing_t m_preview;
  word_parsing_t m_next;

  CifParser(std::istream &io_) : io(io_) { init_fsm(); };

  bool is_newline(char c) { return c == '\r' || c == '\n'; }
  bool is_space(char c) { return c == ' ' || c == '\t'; }
  bool is_period(char c) { return c == '.'; }
  bool is_escape(char c) { return c == '\\'; }
  bool is_digit(char c) { return isdigit(c); }
  bool is_alpha(char c) { return isalpha(c) || c == '_'; }
  bool is_eof(char c) { return c == EOF; }

  bool starts_comment(char c) { return c == '#'; }
  bool starts_category(char c) { return c == '_'; }
  bool starts_attribute(char c) { return c == '.'; }
  bool starts_ss(char c) { return c == '\''; }
  bool starts_ds(char c) { return c == '"'; }
  bool starts_ms(char c) { return c == ';'; }

  bool add_char(fsm_status st) {
    return st == ST_CATEGORY || st == ST_ATTRIBUTE || st == ST_SS || st == ST_DS || st == ST_MS || st == ST_MS_NL ||
           st == ST_VALUE;
  }

  bool end_word(fsm_status st0, fsm_status st1) { return add_char(st0) && !add_char(st1); }

  void init_fsm() {
    fsm[ST_START] = [this](char c) {
      if (is_newline(c)) {
        return ST_START;
      } else if (is_space(c)) {
        return ST_START_;
      } else if (starts_comment(c)) {
        return ST_COMMENT;
      } else if (starts_ms(c)) {
        return ST_MS_START;
      } else if (starts_ss(c)) {
        return ST_SS_START;
      } else if (starts_ds(c)) {
        return ST_DS_START;
      } else if (starts_category(c)) {
        return ST_CATEGORY_START;
      } else
        return ST_VALUE;
    };

    fsm[ST_START_] = [this](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else if (starts_comment(c))
        return ST_COMMENT;
      else if (starts_ms(c))
        return ST_ERROR;
      else if (starts_ss(c))
        return ST_SS_START;
      else if (starts_ds(c))
        return ST_DS_START;
      else if (starts_category(c))
        return ST_CATEGORY_START;
      else
        return ST_VALUE;
    };

    fsm[ST_COMMENT] = [this](char c) {
      if (is_newline(c))
        return ST_START;
      else
        return ST_COMMENT;
    };

    fsm[ST_SS_START] = [this](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else if (is_escape(c))
        return ST_SS_ESCAPE;
      else if (starts_ss(c))
        return ST_ES;
      else
        return ST_SS;
    };

    fsm[ST_SS] = [this](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else if (is_escape(c))
        return ST_SS_ESCAPE;
      else if (starts_ss(c))
        return ST_ES;
      else
        return ST_SS;
    };

    fsm[ST_SS_ESCAPE] = [this](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else
        return ST_SS;
    };

    fsm[ST_DS_START] = [this](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else if (is_escape(c))
        return ST_DS_ESCAPE;
      else if (starts_ds(c))
        return ST_ES;
      else
        return ST_DS;
    };

    fsm[ST_DS] = [this](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else if (is_escape(c))
        return ST_DS_ESCAPE;
      else if (starts_ds(c))
        return ST_ES;
      else
        return ST_DS;
    };

    fsm[ST_DS_ESCAPE] = [this](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else
        return ST_DS;
    };

    fsm[ST_MS_START] = [this](char c) {
      if (is_newline(c))
        return ST_MS_NL;
      else if (is_escape(c))
        return ST_MS_ESCAPE;
      else
        return ST_MS;
    };

    fsm[ST_MS] = [this](char c) {
      if (is_newline(c))
        return ST_MS_NL;
      else if (is_escape(c))
        return ST_MS_ESCAPE;
      else
        return ST_MS;
    };

    fsm[ST_MS_ESCAPE] = [this](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else
        return ST_MS;
    };

    fsm[ST_MS_NL] = [this](char c) {
      if (is_newline(c))
        return ST_MS_NL;
      else if (is_escape(c))
        return ST_MS_ESCAPE;
      else if (starts_ms(c))
        return ST_ES;
      else
        return ST_MS;
    };

    fsm[ST_ES] = [this](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else
        return ST_ERROR;
    };

    fsm[ST_CATEGORY_START] = [this](char c) {
      if (is_newline(c) || is_space(c) || starts_comment(c)) {
        return ST_ERROR;
      } else {
        return ST_CATEGORY;
      }
    };

    fsm[ST_CATEGORY] = [this](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else if (starts_comment(c))
        return ST_COMMENT;
      else if (starts_attribute(c))
        return ST_ATTRIBUTE_START;
      else
        return ST_CATEGORY;
    };

    fsm[ST_ATTRIBUTE_START] = [this](char c) {
      if (is_newline(c) || is_space(c) || starts_comment(c)) {
        return ST_ERROR;
      } else {
        return ST_ATTRIBUTE;
      }
    };

    fsm[ST_ATTRIBUTE] = [this](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else if (starts_comment(c))
        return ST_COMMENT;
      else
        return ST_ATTRIBUTE;
    };

    fsm[ST_VALUE] = [this](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else if (starts_comment(c))
        return ST_COMMENT;
      else
        return ST_VALUE;
    };
  }

  int preview_word() {
    if (!m_preview.empty)
      return 1;

    std::stringstream ss;
    while (true) {
      if (m_preview.status == ST_EOF) {
        return 0;
      } else {
        // std::cout << "status: " << m_preview.status << std::endl;
        m_preview.pre_status = m_preview.status;
        char c = io.get();
        // std::cout << c;
        if (is_eof(c)) {
          m_preview.status = ST_EOF;
        } else {
          m_preview.status = fsm[m_preview.pre_status](c);
          if (m_preview.status == ST_START) {
            m_preview.iline++;
            m_preview.col = 0;
          } else {
            m_preview.col++;
          }
        }
        // std::cout << "status: " << m_preview.status << std::endl;
        if (m_preview.status == ST_ERROR) {
          // std::cout << m_preview.pre_status << ' ' << c << ' ' << m_preview.status << std::endl;
          std::string msg = jnc::string_format("Error: Line %d, Column: %d", m_preview.iline, m_preview.col);
          throw std::runtime_error(msg.c_str());
        } else if (add_char(m_preview.status)) {
          ss << c;
        } else if (add_char(m_preview.pre_status) && !add_char(m_preview.status)) {
          m_preview.word = ss.str();
          if (m_preview.pre_status == ST_CATEGORY)
            m_preview.word_type = WD_CATEGORY;
          else if (m_preview.pre_status == ST_ATTRIBUTE)
            m_preview.word_type = WD_ATTRIBUTE;
          else
            m_preview.word_type = WD_VALUE;
          // std::cout << "Preview: " << m_preview.word << std::endl;
          m_preview.empty = false;
          ss.str("");
          return 1;
        }
      }
    }
  }

  int fetch_word() {
    if (m_preview.empty) {
      if (!preview_word())
        return 0;
    }
    m_next = m_preview;
    m_preview.empty = true;
    return 1;
  }

  void assert_category() {
    if (!fetch_word() || m_next.word_type != WD_CATEGORY) {
      std::string msg = jnc::string_format("Need category near %s: Line %d, Column: %d", m_next.word.c_str(),
                                           m_preview.iline, m_preview.col);
      throw std::runtime_error(msg.c_str());
    }
  }

  void assert_attribute() {
    if (!fetch_word() || m_next.word_type != WD_ATTRIBUTE) {
      std::string msg = jnc::string_format("Need attribute near %s: Line %d, Column: %d", m_next.word.c_str(),
                                           m_preview.iline, m_preview.col);
      throw std::runtime_error(msg.c_str());
    }
  }

  void assert_value() {
    if (!fetch_word() || m_next.word_type != WD_VALUE) {
      std::string msg = jnc::string_format("Need value near %s: Line %d, Column: %d", m_next.word.c_str(),
                                           m_preview.iline, m_preview.col);
      throw std::runtime_error(msg.c_str());
    }
  }

  void read_attributes() {
    m_attributes.clear();
    std::string cat;
    while (true) {
      preview_word();
      if (m_preview.word_type == WD_CATEGORY) {
        assert_category();
        if (cat.empty())
          cat = m_next.word;
        else if (cat != m_next.word)
          throw std::runtime_error("Wrong cif format: incoherent category!");

        assert_attribute();
        m_attributes.push_back(m_next.word);
      } else {
        if (!cat.empty())
          m_category = cat;
        return;
      }
    }
  }

  int next(std::string &cat, std::map<std::string, std::string> &values) {
    if (preview_word()) {
      if (m_preview.word_type == WD_CATEGORY) {
        m_attributes.clear();
        values.clear();

        assert_category();
        cat = m_next.word;
        assert_attribute();
        auto attr = m_next.word;
        assert_value();
        values[attr] = m_next.word;
        return 1;
      } else if (m_preview.word_type == WD_VALUE) {
        if (m_preview.word == "loop_") {
          // std::cout << "LOOP" << std::endl;
          fetch_word();
          read_attributes();
          // std::cout << "Category: " << m_category << std::endl;
          // for (auto && attr : m_attributes) {
          //   std::cout << attr << ' ';
          // }
          // std::cout << std::endl;
          return next(cat, values);
        } else {
          // std::cout << "VALUE" << std::endl;
          if (m_attributes.empty()) {
            while (true) {
              fetch_word();
              preview_word();
              if (m_preview.word_type != WD_VALUE)
                return next(cat, values);
            }
          }
          values.clear();
          cat = m_category;
          for (int i = 0; i < m_attributes.size(); i++) {
            auto &attr = m_attributes[i];
            // std::cout << m_preview.word << ' ' << m_preview.word_type << std::endl;
            assert_value();
            values[attr] = m_next.word;
          }
          return 1;
        }
      }
    } else {
      return 0;
    }
  }
};

} // namespace bio
} // namespace jnc
