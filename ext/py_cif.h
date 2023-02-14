#pragma once

#include "./py_utils.h"
#include <cctype>

namespace jnc {
namespace bio {

class CifParser {
  std::deque<std::string> m_attributes;
  std::string m_category;

  int iline;
  std::istream io;
  int row;
  int col;

  enum FSM_Status {
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
    ST_WORD,
    ST_WORD_PERIOD,
    ST_WORD2,
    ST_PERIOD,
    ST_INT,
    ST_SIGN,
    ST_FRAC_ONLY,
    ST_FRAC,
    ST_SCI_START,
    ST_SCI_SIGN,
    ST_SCI,
    ST_EOF,
    ST_ERROR
  };

  // word parsing results
  struct word_parsing_t {
    std::string word;
    FSM_Status status = ST_START;
    FSM_Status pre_status = ST_START;
  };

  word_parsing_t m_preview;
  word_parsing_t m_next;

  CifParser(std::istream &io_) : io(io_) { init_fsm(); };

  bool is_newline(char c) { return c == '\r' || c == '\n'; }
  bool is_space(char c) { return c == ' ' || c == '\t'; }
  bool is_period(char c) { return c == '.'; }
  bool is_escape(char c) { return '\\'; }
  bool is_digit(char c) { return isdigit(c); }
  bool is_alpha(char c) { return isalpha(c) || c == '_'; }
  bool is_eof(char c) { return c == EOF; }

  bool starts_comment(char c) { return c == '#'; }
  bool starts_category(char c) { return c == '_'; }
  bool starts_attribute(char c) { return c == '.'; }
  bool starts_ss(char c) { return '\''; }
  bool starts_ds(char c) { return '"'; }
  bool starts_ms(char c) { return c == ';'; }
  bool starts_signed(c) { return c == '-' || c == '+'; }
  bool starts_scientific(c) { return c == 'e' || c == 'E'; }

  bool add_char(FSM_Status st) {
    return st == ST_CATEGORY || st == ST_SS || st == ST_DS || st == ST_MS || st == PERIOD || st == ;
  }

  void init_fsm() {
    fsm[ST_START] = [](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else if (starts_comment(c))
        return ST_COMMENT;
      else if (is_period(c))
        return ST_PERIOD;
      else if (starts_ms(c))
        return ST_MS_START;
      else if (starts_ss(c))
        return ST_SS_START;
      else if (starts_ds(c))
        return ST_DS_START;
      else if (starts_category(c))
        return ST_CATEGORY_START;
      else if (is_digit(c))
        return ST_INT;
      else if (is_alpha(c))
        return ST_WORD;
      else if (starts_signed(c))
        return ST_SIGN;
      else if (is_eof(c))
        return ST_EOF;
      else
        return ST_ERROR;
    };

    fsm[ST_START_] = [](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else if (starts_comment(c))
        return ST_COMMENT;
      else if (is_period(c))
        return ST_PERIOD;
      else if (starts_ms(c))
        return ST_ERROR;
      else if (starts_ss(c))
        return ST_SS_START;
      else if (starts_ds(c))
        return ST_DS_START;
      else if (is_digit(c))
        return ST_INT;
      else if (is_alpha(c))
        return ST_WORD;
      else if (starts_signed(c))
        return ST_SIGN;
      else
        return ST_ERROR;
    };

    fsm[ST_COMMENT] = [](char c) {
      if (is_newline(c))
        return ST_START;
      else
        return ST_COMMENT;
    };

    fsm[ST_SS_START] = [](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else if (is_escape(c))
        return ST_SS_ESCAPE;
      else if (starts_ss(c))
        return ST_ES;
      else
        return ST_SS;
    };

    fsm[ST_SS] = [](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else if (is_escape(c))
        return ST_SS_ESCAPE;
      else if (starts_ss(c))
        return ST_ES;
      else
        return ST_SS;
    };

    fsm[ST_SS_ESCAPE] = [](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else
        return ST_SS;
    };

    fsm[ST_DS_START] = [](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else if (is_escape(c))
        return ST_DS_ESCAPE;
      else if (starts_ds(c))
        return ST_ES;
      else
        return ST_DS;
    };

    fsm[ST_DS] = [](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else if (is_escape(c))
        return ST_DS_ESCAPE;
      else if (starts_ds(c))
        return ST_ES;
      else
        return ST_DS;
    };

    fsm[ST_DS_ESCAPE] = [](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else
        return ST_DS;
    };

    fsm[ST_MS_START] = [](char c) {
      if (is_newline(c))
        return ST_MS_NL;
      else if (is_escape(c))
        return ST_MS_ESCAPE;
      else
        return ST_MS;
    };

    fsm[ST_MS] = [](char c) {
      if (is_newline(c))
        return ST_MS_NL;
      else if (is_escape(c))
        return ST_MS_ESCAPE;
      else
        return ST_MS;
    };

    fsm[ST_MS_ESCAPE] = [](char c) {
      if (is_newline(c))
        return ST_ERROR;
      else
        return ST_MS;
    };

    fsm[ST_MS_NL] = [](char c) {
      if (is_newline(c))
        return ST_MS_NL;
      else if (is_escape(c))
        return ST_MS_ESCAPE;
      else if (starts_ms(c))
        return ST_ES;
      else
        return ST_MS;
    };

    fsm[ST_ES] = [](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else
        return ST_ERROR;
    };

    fms[ST_PERIOD] = [](char c) {
      if (is_digit(c))
        return ST_FRAC_ONLY;
      else
        return ST_ERROR;
    };

    fsm[ST_FRAC_ONLY] = [](char c) {
      if (is_digit(c))
        return ST_FRAC_ONLY;
      else if (starts_scientific(c))
        return ST_SCI_START;
      else if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else
        return ST_ERROR;
    };

    fsm[ST_SCI_START] = [](char c) {
      if (starts_signed(c))
        return ST_SCI_SIGN;
      else if (is_digit(c))
        return ST_SCI;
      else
        return ST_ERROR;
    };

    fsm[ST_SCI_SIGN] = [](char c) {
      if (is_digit(c))
        return ST_SCI;
      else
        return ST_ERROR;
    };

    fsm[ST_SCI] = [](char c) {
      if (is_digit(c))
        return ST_SCI;
      else if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else
        return ST_ERROR;
    };

    fsm[ST_INT] = [](char c) {
      if (is_digit(c))
        return ST_INT;
      else if (is_period(c))
        return ST_FRAC;
      else if (starts_scientific(c))
        return ST_SCI_START;
      else if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else
        return ST_ERROR;
    };

    fsm[ST_FRAC] = [](char c) {
      if (is_digit(c))
        return ST_FRAC;
      else if (starts_scientific(c))
        return ST_SCI_START;
      else if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else
        return ST_ERROR;
    };

    fsm[ST_SIGN] = [](char c) {
      if (is_digit(c))
        return ST_INT;
      else if (is_period(c))
        return ST_PERIOD;
      else
        return ST_ERROR;
    };

    fsm[ST_CATEGORY_START] = [](char c) {
      if (is_alpha(c) || is_digit(c))
        return ST_CATEGORY;
      else
        return ST_ERROR;
    };

    fsm[ST_CATEGORY] = [](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else if (is_alpha(c) || is_digit(c))
        return ST_CATEGORY;
      else if (starts_attribute(c))
        return ST_ATTRIBUTE_START;
      else
        return ST_ERROR;
    };

    fsm[ST_ATTRIBUTE_START] = [](char c) {
      if (is_alpha(c) || is_digit(c))
        return ST_ATTRIBUTE;
      else
        return ST_ERROR;
    };

    fsm[ST_ATTRIBUTE] = [](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else if (starts_comment(c))
        return ST_COMMENT;
      else if (is_alpha(c) || is_digit(c))
        return ST_ATTRIBUTE;
      else
        return ST_ERROR;
    };

    fsm[ST_WORD] = [](char c) {
      if (is_newline(c))
        return ST_START;
      else if (is_space(c))
        return ST_START_;
      else if (starts_comment(c))
        return ST_COMMENT;
      else
        return ST_WORD;
    };
  }

  int preview_word() {
    std::stringstream ss;
    while (true) {
      if (m_preview.status == ST_EOF) {
        return 0;
      } else {
        char c = io.get();
        m_preview.pre_status = m_preview.status;
        m_preview.status = fsm[m_preview.status](c);
        if (m_preview.status == ST_ERROR) {
          throw std::runtime_error("Wrong cif format.");
        } else if (add_char(m_preview.status)) {
          ss << c;
        } else if (do_fetch(m_preview.pre_status) && !do_fetch(m_preview.status)) {
          m_preview.word = ss.str();
          ss.str("");
          return 1;
        }
      }
    }
  }

  int fetch_word() {
    if (m_preview.word.empty()) {
      if (!preview())
        return 0;
    }
    m_next = m_preview;
    m_preview.word = "";
    return 1;
  }

  void assert_category() {
    if (!fetch_word() || m_next.pre_status != ST_CATEGORY)
      throw std::runtime_error("Wrong cif format.");
  }

  void assert_attribute() {
    if (!fetch_word() || m_next.pre_status != ST_ATTRIBUTE)
      throw std::runtime_error("Wrong cif format.");
  }

  void assert_word() {
    if (!fetch_word())
      throw std::runtime_error("Wrong cif format.");
  }

  void read_attributes() {
    m_attributes.clear();
    std::string cat;
    while (true) {
      preview_word();
      if (m_preview.pre_status == ST_WORD) {
        assert_category();
        if (cat.empty())
          cat = m_next.word;
        else if (cat != m_next.word)
          throw std::runtime_error("Wrong cif format.");

        assert_attribute();
        m_attributes.push_back(m_next.word)
      } else {
        if (!cat.empty())
          m_category = cat;
        return;
      }
    }
  }

  int next(std::string &cat, std::map<std::string, std::string> &values) {
    if (preview_word()) {
      if (m_preview.pre_status == ST_WORD) {
        if (m_preview.word == "_loop") {
          fetch_word();
          read_attributes();
          return next(cat, values);
        } else {
          values.clear();

          assert_category();
          m_category = m_next.word;
          assert_designate(".");
          assert_attribute();
          auto attr = m_word;

          assert_word();
          values[attr] = m_word;

          return 1;
        }
      } else if (is_value(m_preview.pre_status)) {
        values.clear();
        cat = category;
        values[m_attributes[i]] = m_word;
        for (int i = 1; i < m_attributes.size(); i++) {
          auto &attr = m_attributes[i];
          if (!fetch_word() || !is_value(st0))
            throw std::runtime_error("Wrong cif format.");
          values[attr] = m_word;
        }
        return 1;
      }
    } else {
      return 0;
    }
  }
};

} // namespace bio
} // namespace jnc
