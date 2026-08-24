"""Offline Streamlit smoke tests for the main dashboard workflow."""

from pathlib import Path

from streamlit.testing.v1 import AppTest


APP_PATH = Path(__file__).with_name("app.py")
SAMPLE_SEQUENCE = (
    "ATGGCTAGCTAGCTAGGCTAGCATCGATCGATCGGATCGATCGATCGATCGGCTAGCTAGCTAGCTAGG"
    * 4
)


def test_dashboard_loads_without_exception():
    app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()

    assert not app.exception
    assert any(button.label == "Design and analyze guides" for button in app.button)
    assert any(toggle.label == "Dark mode" for toggle in app.toggle)


def test_pasted_sequence_completes_full_report():
    app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()
    app.radio[0].set_value("Paste sequence").run()
    app.text_area[0].input(SAMPLE_SEQUENCE)
    app.button[0].click().run()

    assert not app.exception
    assert len(app.metric) >= 5
    assert len(app.tabs) == 6
    assert any("Analysis complete" in message.value for message in app.success)


def test_short_sequence_shows_clear_error():
    app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()
    app.radio[0].set_value("Paste sequence").run()
    app.text_area[0].input("ATGC")
    app.button[0].click().run()

    assert not app.exception
    assert any("at least 50 bp" in message.value for message in app.error)

