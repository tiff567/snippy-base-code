# Import smtplib for the actual sending function
import smtplib

# Import the email modules we'll need
from email.message import EmailMessage

def sent_email_notification(from_email_address, to_email_address, subject, msg_filename):

    # Open the plain text file whose name is in textfile for reading.
    with open(msg_filename) as fp:
        # Create a text/plain message
        msg = EmailMessage()
        msg.set_content(fp.read())

    msg['Subject'] = subject
    msg['From'] = from_email_address
    msg['To'] = to_email_address

    # Send the message via our own SMTP server.
    s = smtplib.SMTP('localhost')
    s.send_message(msg)
    s.quit()

    return ("Email notification has been sent to %s", to_email_address)

